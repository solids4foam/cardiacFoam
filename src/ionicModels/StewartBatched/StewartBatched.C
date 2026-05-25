/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "StewartBatched.H"
#include "Stewart_2009Batch.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& StewartFamilyInfo();
}

namespace
{
    enum class StewartRLTauSource : unsigned char
    {
        none,
        lookup
    };

    struct StewartRLDispatchEntry
    {
        StewartRLTauSource tauSource;
        Foam::label tauIndex;
    };

    inline StewartRLDispatchEntry rlNone()
    {
        return {StewartRLTauSource::none, -1};
    }

    inline StewartRLDispatchEntry rlLookup(const Foam::label tauI)
    {
        return {StewartRLTauSource::lookup, tauI};
    }

    const std::array<StewartRLDispatchEntry, NUM_STATES>
        StewartRushLarsenDispatch = []()
    {
        std::array<StewartRLDispatchEntry, NUM_STATES> e{};
        e.fill(rlNone());

        e[Ihyperpolarization_activated_current_y_gate_y] = rlLookup(AV_tau_y);
        e[Irapid_time_dependent_potassium_current_Xr1_gate_Xr1] = rlLookup(AV_tau_xr1);
        e[Irapid_time_dependent_potassium_current_Xr2_gate_Xr2] = rlLookup(AV_tau_xr2);
        e[Islow_time_dependent_potassium_current_Xs_gate_Xs] = rlLookup(AV_tau_xs);
        e[Ifast_sodium_current_m_gate_m] = rlLookup(AV_tau_m);
        e[Ifast_sodium_current_h_gate_h] = rlLookup(AV_tau_h);
        e[Ifast_sodium_current_j_gate_j] = rlLookup(AV_tau_j);
        e[IL_type_Ca_current_d_gate_d] = rlLookup(AV_tau_d);
        e[IL_type_Ca_current_f_gate_f] = rlLookup(AV_tau_f);
        e[IL_type_Ca_current_f2_gate_f2] = rlLookup(AV_tau_f2);
        e[IL_type_Ca_current_fCass_gate_fCass] = rlLookup(AV_tau_fCass);
        e[Itransient_outward_current_s_gate_s] = rlLookup(AV_tau_s);
        e[Itransient_outward_current_r_gate_r] = rlLookup(AV_tau_r);

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
    void launchStewartBatchKernel
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

    void launchStewartEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchStewartScaleIonKernel
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
    defineTypeNameAndDebug(StewartBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        StewartBatched,
        dictionary
    );
}

Foam::StewartBatched::StewartBatched
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
        StewartFamilyInfo()
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
            Info<< "StewartBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "StewartBatched: useSoAEvaluator is on but "
                << "no CUDA device is visible (" << cudaGetErrorString(err)
                << "); falling back to the host SIMD path." << nl;
        }
    }
#endif

    if (useSoAEvaluator_)
    {
        setHotPathSupportSize(NUM_STEWART_BATCH_SUPPORT);

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

    StewartinitConsts
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

Foam::StewartBatched::~StewartBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}




Foam::List<Foam::word> Foam::StewartBatched::supportedTissueTypes() const
{
    return {"myocyte"};
}

void Foam::StewartBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    StewartcomputeVariables
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

void Foam::StewartBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_STEWART_BATCH_SUPPORT),
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


void Foam::StewartBatched::importFields
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


void Foam::StewartBatched::solveODE
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

void Foam::StewartBatched::solveBatched
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
            static_cast<std::size_t>(NUM_STEWART_BATCH_SUPPORT),
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

            launchStewartBatchKernel
            (
                tSub, cuda_.d_constants,
                static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                solveVm, stimulusPOD_
            );

            launchStewartEulerStepKernel
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
        launchStewartBatchKernel
        (
            tEnd, cuda_.d_constants,
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );

        launchStewartScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 1.0,
            static_cast<int>(N),
            static_cast<int>(STEWART_BATCH_SUPPORT_Iion_cm)
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

        StewartComputeVariablesBatch
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
    StewartComputeVariablesBatch
    (
        tEnd, CONSTS, static_cast<int>(N), 0, static_cast<int>(N),
        STATES_SoA, RATES_SoA, SUPPORT_SoA,
        solveVm, stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = STEWART_BATCH_SUPPORT_Iion_cm*N;
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = SUPPORT_SoA[IionBase + cellI];
        }
    }

    setIOEvaluationModelTime(tEnd);
}

bool Foam::StewartBatched::rushLarsenParameters
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

    const StewartRLDispatchEntry& entry = StewartRushLarsenDispatch[stateI];

    switch (entry.tauSource)
    {
        case StewartRLTauSource::lookup:
            tau = algebraicValues[entry.tauIndex];
            break;

        case StewartRLTauSource::none:
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

void Foam::StewartBatched::derivatives
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
