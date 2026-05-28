/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "GrandiBatched.H"
#include "Grandi_2011Batch.H"
#include "batchedRushLarsenEntry.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& GrandiFamilyInfo();
}

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    GrandiRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[Ikr_xr]       = Foam::rlScalarAlgAndSupport(AV_ikr_xr_tau,       AV_ikr_xr_inf,       Foam::GRANDI_BATCH_SUPPORT_tau_Ikr_xr,       Foam::GRANDI_BATCH_SUPPORT_gInf_Ikr_xr);
        t[Iks_xs]       = Foam::rlScalarAlgAndSupport(AV_iks_xs_tau,        AV_iks_xs_inf,        Foam::GRANDI_BATCH_SUPPORT_tau_Iks_xs,       Foam::GRANDI_BATCH_SUPPORT_gInf_Iks_xs);
        t[Ikur_ikur_r]  = Foam::rlScalarAlgAndSupport(AV_ikur_r_tau,        AV_ikur_r_inf,        Foam::GRANDI_BATCH_SUPPORT_tau_Ikur_ikur_r,  Foam::GRANDI_BATCH_SUPPORT_gInf_Ikur_ikur_r);
        t[Ikur_s]       = Foam::rlScalarAlgAndSupport(AV_ikur_s_tau,        AV_ikur_s_inf,        Foam::GRANDI_BATCH_SUPPORT_tau_Ikur_s,       Foam::GRANDI_BATCH_SUPPORT_gInf_Ikur_s);
        t[Ina_h]        = Foam::rlScalarAlgAndSupport(AV_ina_h_tau,         AV_ina_h_inf,         Foam::GRANDI_BATCH_SUPPORT_tau_Ina_h,        Foam::GRANDI_BATCH_SUPPORT_gInf_Ina_h);
        t[Ina_j]        = Foam::rlScalarAlgAndSupport(AV_ina_j_tau,         AV_ina_j_inf,         Foam::GRANDI_BATCH_SUPPORT_tau_Ina_j,        Foam::GRANDI_BATCH_SUPPORT_gInf_Ina_j);
        t[Ina_m]        = Foam::rlScalarAlgAndSupport(AV_ina_m_tau,         AV_ina_m_inf,         Foam::GRANDI_BATCH_SUPPORT_tau_Ina_m,        Foam::GRANDI_BATCH_SUPPORT_gInf_Ina_m);
        t[Inal_hl]      = Foam::rlScalarConstTauAlgInfAndSupport(AC_inal_hl_tau, AV_inal_hl_inf, Foam::GRANDI_BATCH_SUPPORT_tau_Inal_hl, Foam::GRANDI_BATCH_SUPPORT_gInf_Inal_hl);
        t[Inal_ml]      = Foam::rlScalarAlgAndSupport(AV_inal_ml_tau,       AV_inal_ml_inf,       Foam::GRANDI_BATCH_SUPPORT_tau_Inal_ml,      Foam::GRANDI_BATCH_SUPPORT_gInf_Inal_ml);
        t[Ical_d]       = Foam::rlScalarAlgAndSupport(AV_ical_d_tau,        AV_ical_d_inf,        Foam::GRANDI_BATCH_SUPPORT_tau_Ical_d,       Foam::GRANDI_BATCH_SUPPORT_gInf_Ical_d);
        t[Ical_f]       = Foam::rlScalarAlgAndSupport(AV_ical_f_tau,        AV_ical_f_inf,        Foam::GRANDI_BATCH_SUPPORT_tau_Ical_f,       Foam::GRANDI_BATCH_SUPPORT_gInf_Ical_f);
        t[Ical_fCaB_jn] = Foam::rlScalarAlgAndSupport(AV_ical_fCaB_jn_tau, AV_ical_fCaB_jn_inf, Foam::GRANDI_BATCH_SUPPORT_tau_Ical_fCaB_jn, Foam::GRANDI_BATCH_SUPPORT_gInf_Ical_fCaB_jn);
        t[Ical_fCaB_sl] = Foam::rlScalarAlgAndSupport(AV_ical_fCaB_sl_tau, AV_ical_fCaB_sl_inf, Foam::GRANDI_BATCH_SUPPORT_tau_Ical_fCaB_sl, Foam::GRANDI_BATCH_SUPPORT_gInf_Ical_fCaB_sl);
        t[Ito_x]        = Foam::rlScalarAlgAndSupport(AV_ito_x_tau,         AV_ito_x_inf,         Foam::GRANDI_BATCH_SUPPORT_tau_Ito_x,        Foam::GRANDI_BATCH_SUPPORT_gInf_Ito_x);
        t[Ito_y]        = Foam::rlScalarAlgAndSupport(AV_ito_y_tau,         AV_ito_y_inf,         Foam::GRANDI_BATCH_SUPPORT_tau_Ito_y,        Foam::GRANDI_BATCH_SUPPORT_gInf_Ito_y);

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
    void launchGrandiBatchKernel
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

    void launchGrandiEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchGrandiScaleIonKernel
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
    bool useGrandiCompactSupport(const dictionary& dict)
    {
        const word modelName =
            dict.lookupOrDefault<word>("ionicModel", word::null);

        return modelName == "GrandicompactBatched"
            || dict.lookupOrDefault<Switch>("useCompactSupport", false);
    }

    defineTypeNameAndDebug(GrandiBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        GrandiBatched,
        dictionary
    );
    defineTypeNameAndDebug(GrandicompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        GrandicompactBatched,
        dictionary
    );
}

Foam::GrandiBatched::GrandiBatched
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
        GrandiFamilyInfo()
    ),
    useSoAEvaluator_
    (
        dict.lookupOrDefault<Switch>("useSoAEvaluator", false)
    ),
    useCompactSupport_(useGrandiCompactSupport(dict)),
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
            Info<< "GrandiBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "GrandiBatched: useSoAEvaluator is on but "
                << "no CUDA device is visible (" << cudaGetErrorString(err)
                << "); falling back to the host SIMD path." << nl;
        }
    }
#endif

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    GrandiinitConsts
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

    if (useSoAEvaluator_ || useCompactSupport_)
    {
        setHotPathSupportSize(NUM_GRANDI_BATCH_SUPPORT);
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
                << "` was requested. The SoA path currently uses "
                << "explicit Euler. Set `batchedIntegrator euler;` "
                << "to silence this warning, or unset "
                << "`useSoAEvaluator` to use the cell-major path."
                << nl;
        }
    }
}

Foam::GrandicompactBatched::GrandicompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    GrandiBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

Foam::GrandiBatched::~GrandiBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}




Foam::List<Foam::word> Foam::GrandiBatched::supportedTissueTypes() const
{
    return {"myocyte"};
}

void Foam::GrandiBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    GrandicomputeVariables
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

Foam::scalar Foam::GrandiBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[GRANDI_BATCH_SUPPORT_Iion_cm];
}

void Foam::GrandiBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (const auto& entry : GrandiRushLarsenDispatch)
    {
        projectScalarRushLarsenEntryToSupport
        (
            entry,
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }

    supportValues[GRANDI_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}

bool Foam::GrandiBatched::rushLarsenParametersFromHotPathSupport
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

    const auto& entry = GrandiRushLarsenDispatch[stateI];
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

void Foam::GrandiBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_GRANDI_BATCH_SUPPORT),
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


void Foam::GrandiBatched::importFields
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


void Foam::GrandiBatched::solveODE
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

void Foam::GrandiBatched::solveBatched
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
            static_cast<std::size_t>(NUM_GRANDI_BATCH_SUPPORT),
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

            launchGrandiBatchKernel
            (
                tSub, cuda_.d_constants,
                static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                solveVm, stimulusPOD_
            );

            launchGrandiEulerStepKernel
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
        launchGrandiBatchKernel
        (
            tEnd, cuda_.d_constants,
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );

        launchGrandiScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 1.0,
            static_cast<int>(N),
            static_cast<int>(GRANDI_BATCH_SUPPORT_Iion_cm)
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

        GrandiComputeVariablesBatch
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
    GrandiComputeVariablesBatch
    (
        tEnd, CONSTS, static_cast<int>(N), 0, static_cast<int>(N),
        STATES_SoA, RATES_SoA, SUPPORT_SoA,
        solveVm, stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = GRANDI_BATCH_SUPPORT_Iion_cm*N;
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = SUPPORT_SoA[IionBase + cellI];
        }
    }

    setIOEvaluationModelTime(tEnd);
}

bool Foam::GrandiBatched::rushLarsenParameters
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

    const auto& entry = GrandiRushLarsenDispatch[stateI];
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

void Foam::GrandiBatched::derivatives
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
