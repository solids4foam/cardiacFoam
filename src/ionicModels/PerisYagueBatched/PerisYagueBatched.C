/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "PerisYagueBatched.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& PerisYagueFamilyInfo();
}

namespace
{
    enum class PerisYagueRLTauSource : unsigned char
    {
        none,
        lookup,
        constant,
        literal
    };

    struct PerisYagueRLDispatchEntry
    {
        PerisYagueRLTauSource tauSource;
        Foam::label tauIndex;
        Foam::scalar tauLiteral;
    };

    inline PerisYagueRLDispatchEntry rlNone()
    {
        return {PerisYagueRLTauSource::none, -1, 0.0};
    }

    inline PerisYagueRLDispatchEntry rlLookup(const Foam::label tauI)
    {
        return {PerisYagueRLTauSource::lookup, tauI, 0.0};
    }

    inline PerisYagueRLDispatchEntry rlConstant(const Foam::label tauI)
    {
        return {PerisYagueRLTauSource::constant, tauI, 0.0};
    }

    inline PerisYagueRLDispatchEntry rlLiteral(const Foam::scalar tau)
    {
        return {PerisYagueRLTauSource::literal, -1, tau};
    }

    const std::array<PerisYagueRLDispatchEntry, NUM_STATES>
        PerisYagueRushLarsenDispatch = []()
    {
        std::array<PerisYagueRLDispatchEntry, NUM_STATES> e{};
        e.fill(rlNone());

        e[ina_m]     = rlLookup(AV_ina_m_tau);
        e[ina_h]     = rlLookup(AV_ina_h_tau);
        e[ina_j]     = rlLookup(AV_ina_j_tau);
        e[ikr_xr]    = rlLookup(AV_ikr_xr_tau);
        e[iks_xs]    = rlLookup(AV_iks_xs_tau);
        e[ikur_ua]   = rlLookup(AV_ikur_ua_tau);
        e[ikur_uif]  = rlLookup(AV_ikur_uif_tau);
        e[ikur_uis]  = rlLookup(AV_ikur_uis_tau);
        e[ical_d]    = rlLookup(AV_ical_d_tau);
        e[ical_f]    = rlLookup(AV_ical_f_tau);
        e[ical_fCa]  = rlConstant(AC_ical_fCa_tau);
        e[iclca_qCa] = rlLiteral(2.0);
        e[ryr_u]     = rlConstant(AC_cajsr_u_tau);
        e[ryr_w]     = rlLookup(AV_ryr_w_tau);

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
    void launchPerisYagueBatchKernel
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

    void launchPerisYagueEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchPerisYagueScaleIonKernel
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
    defineTypeNameAndDebug(PerisYagueBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        PerisYagueBatched,
        dictionary
    );
}

Foam::PerisYagueBatched::PerisYagueBatched
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
        PerisYagueFamilyInfo()
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
            Info<< "PerisYagueBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "PerisYagueBatched: useSoAEvaluator is on but "
                << "no CUDA device is visible (" << cudaGetErrorString(err)
                << "); falling back to the host SIMD path." << nl;
        }
    }
#endif

    if (useSoAEvaluator_)
    {
        setHotPathSupportSize(NUM_PERISYAGUE_BATCH_SUPPORT);

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

    PerisYague_2022initConsts
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

Foam::PerisYagueBatched::~PerisYagueBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}




Foam::List<Foam::word> Foam::PerisYagueBatched::supportedTissueTypes() const
{
    return {"myocyte"};
}

void Foam::PerisYagueBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    PerisYague_2022computeVariables
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

void Foam::PerisYagueBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_PERISYAGUE_BATCH_SUPPORT),
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


void Foam::PerisYagueBatched::importFields
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


void Foam::PerisYagueBatched::solveODE
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

void Foam::PerisYagueBatched::solveBatched
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    const label N = nCells();
    const label nSub = nSubsteps();
    const scalar dtModel = deltaT * timeScaleFactor();
    const scalar dtSubstep = dtModel / scalar(nSub);
    const scalar tStart = stepStartTime * timeScaleFactor();
    const bool   solveVm = solveVmWithinODESolver();

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

    scalar* const STATES_SoA  = statesSoAData();
    scalar* const RATES_SoA   = ratesSoAData();
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
            static_cast<std::size_t>(NUM_PERISYAGUE_BATCH_SUPPORT),
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

            launchPerisYagueBatchKernel
            (
                tSub, cuda_.d_constants,
                static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                solveVm, stimulusPOD_
            );

            launchPerisYagueEulerStepKernel
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
        launchPerisYagueBatchKernel
        (
            tEnd, cuda_.d_constants,
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );

        launchPerisYagueScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 1.0,
            static_cast<int>(N),
            static_cast<int>(PERISYAGUE_BATCH_SUPPORT_Iion_cm)
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

        PerisYagueComputeVariablesBatch
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
                    dtSubstep * RATES_SoA[base + cellI];
            }
        }
    }

    // Final evaluation at end-of-step to refresh RATES/SUPPORT for the
    // downstream coupling step (Im extraction below).
    const scalar tEnd = tStart + dtModel;
    PerisYagueComputeVariablesBatch
    (
        tEnd, CONSTS,
        static_cast<int>(N), 0, static_cast<int>(N),
        STATES_SoA, RATES_SoA, SUPPORT_SoA,
        solveVm, stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = PERISYAGUE_BATCH_SUPPORT_Iion_cm * N;
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

bool Foam::PerisYagueBatched::rushLarsenParameters
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

    if (stateI == ryr_v)
    {
        const scalar Irel =
            CONSTANTS_[AC_krel]
           *stateValues[ryr_u]*stateValues[ryr_u]
           *stateValues[ryr_v]
           *stateValues[ryr_w]
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

    const PerisYagueRLDispatchEntry& entry =
        PerisYagueRushLarsenDispatch[stateI];

    switch (entry.tauSource)
    {
        case PerisYagueRLTauSource::lookup:
            tau = algebraicValues[entry.tauIndex];
            break;

        case PerisYagueRLTauSource::constant:
            tau = CONSTANTS_[entry.tauIndex];
            break;

        case PerisYagueRLTauSource::literal:
            tau = entry.tauLiteral;
            break;

        case PerisYagueRLTauSource::none:
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

void Foam::PerisYagueBatched::derivatives
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
