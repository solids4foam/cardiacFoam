/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "CourtemancheBatched.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& CourtemancheFamilyInfo();
}

namespace
{
    enum class CourtemancheRLTauSource : unsigned char
    {
        none,
        lookup,
        constant
    };

    struct CourtemancheRLDispatchEntry
    {
        CourtemancheRLTauSource tauSource;
        Foam::label tauIndex;
    };

    inline CourtemancheRLDispatchEntry rlNone()
    {
        return {CourtemancheRLTauSource::none, -1};
    }

    inline CourtemancheRLDispatchEntry rlLookup(const Foam::label tauI)
    {
        return {CourtemancheRLTauSource::lookup, tauI};
    }

    inline CourtemancheRLDispatchEntry rlConstant(const Foam::label tauI)
    {
        return {CourtemancheRLTauSource::constant, tauI};
    }

    const std::array<CourtemancheRLDispatchEntry, NUM_STATES>
        CourtemancheRushLarsenDispatch = []()
    {
        std::array<CourtemancheRLDispatchEntry, NUM_STATES> e{};
        e.fill(rlNone());

        e[ina_m]   = rlLookup(AV_ina_m_tau);
        e[ina_h]   = rlLookup(AV_ina_h_tau);
        e[ina_j]   = rlLookup(AV_ina_j_tau);
        e[ical_d]  = rlLookup(AV_ical_d_tau);
        e[ical_f]  = rlLookup(AV_ical_f_tau);
        e[ical_fCa] = rlConstant(AC_ical_fCa_tau);
        e[ito_oa]  = rlLookup(AV_ito_oa_tau);
        e[ito_oi]  = rlLookup(AV_ito_oi_tau);
        e[ikur_ua] = rlLookup(AV_ikur_ua_tau);
        e[ikur_ui] = rlLookup(AV_ikur_ui_tau);
        e[ikr_xr]  = rlLookup(AV_ikr_xr_tau);
        e[iks_xs]  = rlLookup(AV_iks_xs_tau);
        e[cajsr_u] = rlConstant(AC_cajsr_u_tau);
        e[cajsr_w] = rlLookup(AV_cajsr_w_tau);

        return e;
    }();
}

#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "stimulusIO.H"

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchCourtemancheBatchKernel
    (
        double t,
        const double* d_CONSTANTS,
        int N,
        const double* d_STATES,
        double* d_RATES,
        double* d_SUPPORT,
        bool solveVm,
        StimulusProtocolPOD stimulus
    );

    void launchCourtemancheEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchCourtemancheScaleIonKernel
    (
        const double* d_SUPPORT,
        double* d_Im,
        double scale,
        int N,
        int IionSlot
    );
}

#endif // HAS_CUDA

namespace Foam
{
    defineTypeNameAndDebug(CourtemancheBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        CourtemancheBatched,
        dictionary
    );
}

Foam::CourtemancheBatched::CourtemancheBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    configuredBatchedIonicModel
    (
        dict,
        num,
        initialDeltaT,
        solveVmWithinODESolver,
        NUM_STATES,
        NUM_ALGEBRAIC,
        CourtemancheFamilyInfo()
    ),
    useSoAEvaluator_
    (
        dict.lookupOrDefault<Switch>("useSoAEvaluator", false)
    ),
#ifdef HAS_CUDA
    useDevice_(false),
#endif
    stimulusPOD_()
{
    ionicModel::setTissueFromDict();

#ifdef HAS_CUDA
    if (useSoAEvaluator_)
    {
        int nDevices = 0;
        cudaError_t err = cudaGetDeviceCount(&nDevices);
        if (err == cudaSuccess && nDevices > 0)
        {
            const int rank = Pstream::myProcNo();
            const int devId = rank % nDevices;
            CARDIAC_CUDA_CHECK(cudaSetDevice(devId));
            useDevice_ = true;
            Info<< "CourtemancheBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "CourtemancheBatched: useSoAEvaluator is on but "
                << "no CUDA device is visible (" << cudaGetErrorString(err)
                << "); falling back to the host SIMD path." << nl;
        }
    }
#endif

    if (useSoAEvaluator_)
    {
        setHotPathSupportSize(NUM_COURTEMANCHE_BATCH_SUPPORT);

        const word integrator =
            dict.lookupOrDefault<word>("batchedIntegrator", "euler");
        if (integrator != "euler")
        {
            WarningInFunction
                << "useSoAEvaluator is enabled for "
                << type() << " but `batchedIntegrator " << integrator
                << "` is not yet honored on the batched path; reverting "
                << "to explicit Euler for the SoA solver. Set "
                << "`batchedIntegrator euler;` to silence this warning, "
                << "or unset `useSoAEvaluator` to use the cell-major "
                << "integrator path." << nl;
        }
    }

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    CourtemancheinitConsts
    (
        CONSTANTS_.data(),
        initialRates,
        initialStates,
        tissue(),
        dict
    );

    for (label cellI = 0; cellI < nCells(); ++cellI)
    {
        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            state(cellI, stateI) = initialStates[stateI];
            rate(cellI, stateI) = initialRates[stateI];
        }
    }

    configurePersistentAlgebraics();
    syncAllToIO();

    if (!utilitiesMode())
    {
        setStimulusProtocolFromDict(dict);
    }
}

Foam::CourtemancheBatched::~CourtemancheBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}



Foam::List<Foam::word> Foam::CourtemancheBatched::supportedTissueTypes() const
{
    return {"myocyte"};
}

void Foam::CourtemancheBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    CourtemanchecomputeVariables
    (
        modelTime,
        CONSTANTS_.data(),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        algebraicValues.data(),
        tissue(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

void Foam::CourtemancheBatched::prepareIOAccess
(
    const wordList& requestedNames,
    const bool needsAlgebraics
) const
{
#ifdef HAS_CUDA
    if (useDevice_ && cuda_.deviceDirty)
    {
        cuda_.syncStatesDeviceToHost
        (
            statesSoAData(),
            static_cast<std::size_t>(NUM_STATES),
            static_cast<std::size_t>(nCells())
        );
        cuda_.syncSupportDeviceToHost
        (
            supportSoAData(),
            static_cast<std::size_t>(NUM_COURTEMANCHE_BATCH_SUPPORT),
            static_cast<std::size_t>(nCells())
        );
    }
#endif

    configuredBatchedIonicModel::prepareIOAccess
    (
        requestedNames,
        needsAlgebraics
    );
}


void Foam::CourtemancheBatched::importFields
(
    const volScalarField& Vm,
    const wordList& fieldNames,
    const PtrList<volScalarField>& inFields
)
{
    configuredBatchedIonicModel::importFields(Vm, fieldNames, inFields);

#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.hostDirty = true;
    }
#endif
}


void Foam::CourtemancheBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    if (useSoAEvaluator_ && !utilitiesMode())
    {
        solveBatched(stepStartTime, deltaT, Vm, Im);
        return;
    }

    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

void Foam::CourtemancheBatched::solveBatched
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    const label N = nCells();
    const label nSub = nSubsteps();
    const scalar dtModel = deltaT*timeScaleFactor();
    const scalar dtSubstep = dtModel/scalar(nSub);
    const scalar tStart = stepStartTime*timeScaleFactor();
    const bool solveVm = solveVmWithinODESolver();

    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());

    if (!solveVm)
    {
        for (label cellI = 0; cellI < N; ++cellI)
        {
            state(cellI, membrane_V) = vmToState(Vm[cellI]);
        }
#ifdef HAS_CUDA
        if (useDevice_)
        {
            cuda_.hostDirty = true;
        }
#endif
    }

    scalar* const STATES_SoA = statesSoAData();
    scalar* const RATES_SoA = ratesSoAData();
    scalar* const SUPPORT_SoA = supportSoAData();
    const scalar* const CONSTS = CONSTANTS_.cdata();

    markIODirty();
    setIOEvaluationModelTime(tStart);

#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.allocate
        (
            static_cast<std::size_t>(N),
            static_cast<std::size_t>(NUM_STATES),
            static_cast<std::size_t>(NUM_COURTEMANCHE_BATCH_SUPPORT),
            static_cast<std::size_t>(CONSTANTS_.size())
        );
        cuda_.uploadConstants(CONSTANTS_.cdata(), CONSTANTS_.size());

        if (cuda_.hostDirty)
        {
            cuda_.syncStatesHostToDevice
            (
                statesSoAData(),
                static_cast<std::size_t>(NUM_STATES),
                static_cast<std::size_t>(N)
            );
        }

        for (label sub = 0; sub < nSub; ++sub)
        {
            const scalar tSub = tStart + scalar(sub)*dtSubstep;

            launchCourtemancheBatchKernel
            (
                tSub, cuda_.d_constants,
                static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                solveVm, stimulusPOD_
            );

            launchCourtemancheEulerStepKernel
            (
                cuda_.d_states, cuda_.d_rates,
                static_cast<double>(dtSubstep),
                static_cast<int>(N),
                static_cast<int>(NUM_STATES),
                solveVm,
                static_cast<int>(membrane_V)
            );
        }

        const scalar tEnd = tStart + dtModel;
        launchCourtemancheBatchKernel
        (
            tEnd, cuda_.d_constants,
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );

        launchCourtemancheScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 1.0,
            static_cast<int>(N),
            static_cast<int>(COURTEMANCHE_BATCH_SUPPORT_Iion_cm)
        );

        cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));

        cuda_.deviceDirty = true;
        setIOEvaluationModelTime(tEnd);
        return;
    }
#endif // HAS_CUDA

    for (label sub = 0; sub < nSub; ++sub)
    {
        const scalar tSub = tStart + scalar(sub)*dtSubstep;

        CourtemancheComputeVariablesBatch
        (
            tSub, CONSTS, static_cast<int>(N), 0, static_cast<int>(N),
            STATES_SoA, RATES_SoA, SUPPORT_SoA,
            solveVm, stimulusPOD_
        );

        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            if (!solveVm && stateI == membrane_V) continue;
            const label base = stateI*N;
#ifdef _OPENMP
            #pragma omp simd
#endif
            for (label cellI = 0; cellI < N; ++cellI)
            {
                STATES_SoA[base + cellI] +=
                    dtSubstep*RATES_SoA[base + cellI];
            }
        }
    }

    const scalar tEnd = tStart + dtModel;
    CourtemancheComputeVariablesBatch
    (
        tEnd, CONSTS, static_cast<int>(N), 0, static_cast<int>(N),
        STATES_SoA, RATES_SoA, SUPPORT_SoA,
        solveVm, stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = COURTEMANCHE_BATCH_SUPPORT_Iion_cm*N;
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = SUPPORT_SoA[IionBase + cellI];
        }
    }
    else
    {
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = 0.0;
        }
    }

    setIOEvaluationModelTime(tEnd);
}

bool Foam::CourtemancheBatched::rushLarsenParameters
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& algebraicValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI < 0 || stateI >= NUM_STATES)
    {
        return false;
    }

    if (stateI == cajsr_v)
    {
        const scalar Irel =
            CONSTANTS_[AC_K_rel]
           *stateValues[cajsr_u]*stateValues[cajsr_u]
           *stateValues[cajsr_v]
           *stateValues[cajsr_w]
           *(stateValues[calcium_CaRel] - stateValues[calcium_Cai]);

        const scalar Fn =
            1e-12*CONSTANTS_[AC_V_rel]*Irel
          - 5e-13/CONSTANTS_[AC_F]
           *(0.5*algebraicValues[AV_ICaL] - 0.2*algebraicValues[AV_INaCa])
           *CONSTANTS_[AC_Cm];

        steadyState =
            1.0 - 1.0
           /(1.0 + std::exp(-(Fn - 0.2*CONSTANTS_[AC_c1])/CONSTANTS_[AC_c2]));

        tau =
            1.91 + 2.09
           /(1.0 + std::exp(-(Fn - CONSTANTS_[AC_c1])/CONSTANTS_[AC_c2]));

        return tau > VSMALL
            && std::isfinite(tau)
            && std::isfinite(steadyState);
    }

    const CourtemancheRLDispatchEntry& entry =
        CourtemancheRushLarsenDispatch[stateI];

    switch (entry.tauSource)
    {
        case CourtemancheRLTauSource::lookup:
            tau = algebraicValues[entry.tauIndex];
            break;

        case CourtemancheRLTauSource::constant:
            tau = CONSTANTS_[entry.tauIndex];
            break;

        case CourtemancheRLTauSource::none:
        default:
            return false;
    }

    if (tau <= VSMALL)
    {
        return false;
    }

    steadyState = stateValues[stateI] + rateValues[stateI]*tau;
    return std::isfinite(steadyState) && std::isfinite(tau);
}

void Foam::CourtemancheBatched::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(t, y, dydt, algebraics);
}

// ************************************************************************* //
