/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "StewartBatched.H"
#include "Stewart_2009Batch.H"
#include "batchedRushLarsenEntry.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& StewartFamilyInfo();
}

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    StewartRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[Ihyperpolarization_activated_current_y_gate_y]          = Foam::rlScalarAlgAndSupport(AV_tau_y,     AV_y_inf,     Foam::STEWART_BATCH_SUPPORT_tau_y,     Foam::STEWART_BATCH_SUPPORT_gInf_y);
        t[Irapid_time_dependent_potassium_current_Xr1_gate_Xr1]   = Foam::rlScalarAlgAndSupport(AV_tau_xr1,   AV_xr1_inf,   Foam::STEWART_BATCH_SUPPORT_tau_xr1,   Foam::STEWART_BATCH_SUPPORT_gInf_xr1);
        t[Irapid_time_dependent_potassium_current_Xr2_gate_Xr2]   = Foam::rlScalarAlgAndSupport(AV_tau_xr2,   AV_xr2_inf,   Foam::STEWART_BATCH_SUPPORT_tau_xr2,   Foam::STEWART_BATCH_SUPPORT_gInf_xr2);
        t[Islow_time_dependent_potassium_current_Xs_gate_Xs]      = Foam::rlScalarAlgAndSupport(AV_tau_xs,    AV_xs_inf,    Foam::STEWART_BATCH_SUPPORT_tau_xs,    Foam::STEWART_BATCH_SUPPORT_gInf_xs);
        t[Ifast_sodium_current_m_gate_m]                          = Foam::rlScalarAlgAndSupport(AV_tau_m,     AV_m_inf,     Foam::STEWART_BATCH_SUPPORT_tau_m,     Foam::STEWART_BATCH_SUPPORT_gInf_m);
        t[Ifast_sodium_current_h_gate_h]                          = Foam::rlScalarAlgAndSupport(AV_tau_h,     AV_h_inf,     Foam::STEWART_BATCH_SUPPORT_tau_h,     Foam::STEWART_BATCH_SUPPORT_gInf_h);
        t[Ifast_sodium_current_j_gate_j]                          = Foam::rlScalarAlgAndSupport(AV_tau_j,     AV_j_inf,     Foam::STEWART_BATCH_SUPPORT_tau_j,     Foam::STEWART_BATCH_SUPPORT_gInf_j);
        t[IL_type_Ca_current_d_gate_d]                            = Foam::rlScalarAlgAndSupport(AV_tau_d,     AV_d_inf,     Foam::STEWART_BATCH_SUPPORT_tau_d,     Foam::STEWART_BATCH_SUPPORT_gInf_d);
        t[IL_type_Ca_current_f_gate_f]                            = Foam::rlScalarAlgAndSupport(AV_tau_f,     AV_f_inf,     Foam::STEWART_BATCH_SUPPORT_tau_f,     Foam::STEWART_BATCH_SUPPORT_gInf_f);
        t[IL_type_Ca_current_f2_gate_f2]                          = Foam::rlScalarAlgAndSupport(AV_tau_f2,    AV_f2_inf,    Foam::STEWART_BATCH_SUPPORT_tau_f2,    Foam::STEWART_BATCH_SUPPORT_gInf_f2);
        t[IL_type_Ca_current_fCass_gate_fCass]                    = Foam::rlScalarAlgAndSupport(AV_tau_fCass, AV_fCass_inf, Foam::STEWART_BATCH_SUPPORT_tau_fCass, Foam::STEWART_BATCH_SUPPORT_gInf_fCass);
        t[Itransient_outward_current_s_gate_s]                    = Foam::rlScalarAlgAndSupport(AV_tau_s,     AV_s_inf,     Foam::STEWART_BATCH_SUPPORT_tau_s,     Foam::STEWART_BATCH_SUPPORT_gInf_s);
        t[Itransient_outward_current_r_gate_r]                    = Foam::rlScalarAlgAndSupport(AV_tau_r,     AV_r_inf,     Foam::STEWART_BATCH_SUPPORT_tau_r,     Foam::STEWART_BATCH_SUPPORT_gInf_r);

        return t;
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
    bool useStewartCompactSupport(const dictionary& dict)
    {
        const word modelName =
            dict.lookupOrDefault<word>("ionicModel", word::null);

        return modelName == "StewartcompactBatched"
            || dict.lookupOrDefault<Switch>("useCompactSupport", false);
    }

    defineTypeNameAndDebug(StewartBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        StewartBatched,
        dictionary
    );
    defineTypeNameAndDebug(StewartcompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        StewartcompactBatched,
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
    useCompactSupport_(useStewartCompactSupport(dict)),
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

    if (useSoAEvaluator_ || useCompactSupport_)
    {
        setHotPathSupportSize(NUM_STEWART_BATCH_SUPPORT);
    }

    if (useSoAEvaluator_)
    {
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

Foam::StewartcompactBatched::StewartcompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    StewartBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

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

Foam::scalar Foam::StewartBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[STEWART_BATCH_SUPPORT_Iion_cm];
}

void Foam::StewartBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (const auto& entry : StewartRushLarsenDispatch)
    {
        projectScalarRushLarsenEntryToSupport
        (
            entry,
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }

    supportValues[STEWART_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}

bool Foam::StewartBatched::rushLarsenParametersFromHotPathSupport
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& supportValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI < 0 || stateI >= NUM_STATES) return false;

    const auto& entry = StewartRushLarsenDispatch[stateI];
    return resolveSupportRushLarsenEntry
    (
        entry,
        CONSTANTS_,
        supportValues,
        VSMALL,
        steadyState,
        tau
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
    if (stateI < 0 || stateI >= NUM_STATES) return false;

    const auto& entry = StewartRushLarsenDispatch[stateI];
    return resolveScalarRushLarsenEntry
    (
        entry,
        CONSTANTS_,
        algebraicValues,
        VSMALL,
        steadyState,
        tau
    );
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
