/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
\*---------------------------------------------------------------------------*/

#include "GaurBatched.H"
#include "Gaur_2021Batch.H"
#include "Gaur_2021Names.H"
#include "batchedRushLarsenEntry.H"

namespace Foam
{
    const ionicModelFamilyInfo& GaurFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            GaurCONSTANTS_NAMES,
            GaurSTATES_NAMES,
            GaurALGEBRAIC_NAMES,
            cell_v,
            1000.0,
            1000.0,
            0.0,
            nullptr
        };
        return info;
    }
}
#include "addToRunTimeSelectionTable.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "Pstream.H"
#include <array>

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchGaurBatchKernel
    (
        double t, const double* d_CONSTANTS, int N,
        const double* d_STATES, double* d_RATES, double* d_SUPPORT,
        int tissueFlag, bool solveVm, StimulusProtocolPOD stimulus
    );
    void launchGaurEulerStepKernel
    (
        double* d_STATES, const double* d_RATES, double dt,
        int N, int nStates, bool solveVm, int vmStateI
    );
    void launchGaurScaleIonKernel
    (
        const double* d_SUPPORT, double* d_Im,
        double scale, int N, int IionBase
    );
}

#endif // HAS_CUDA

namespace Foam
{
    bool useGaurCompactSupport(const dictionary& dict)
    {
        const word modelName =
            dict.lookupOrDefault<word>("ionicModel", word::null);

        return modelName == "GaurcompactBatched"
            || dict.lookupOrDefault<Switch>("useCompactSupport", false);
    }

    defineTypeNameAndDebug(GaurBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        GaurBatched,
        dictionary
    );
    defineTypeNameAndDebug(GaurcompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        GaurcompactBatched,
        dictionary
    );
}


namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    GaurRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[I_Na_m]    = Foam::rlScalarAlgAndSupport(AV_tau_m,    AV_m_inf,    Foam::GAUR_BATCH_SUPPORT_tau_m,    Foam::GAUR_BATCH_SUPPORT_gInf_m);
        t[I_Na_h]    = Foam::rlScalarAlgAndSupport(AV_tau_h,    AV_h_inf,    Foam::GAUR_BATCH_SUPPORT_tau_h,    Foam::GAUR_BATCH_SUPPORT_gInf_h);
        t[I_Na_j]    = Foam::rlScalarAlgAndSupport(AV_tau_j,    AV_j_inf,    Foam::GAUR_BATCH_SUPPORT_tau_j,    Foam::GAUR_BATCH_SUPPORT_gInf_j);

        t[INaL_ml]   = Foam::rlScalarAlgAndSupport(AV_tau_ml,   AV_ml_inf,   Foam::GAUR_BATCH_SUPPORT_tau_ml,   Foam::GAUR_BATCH_SUPPORT_gInf_ml);
        t[INaL_hl]   = Foam::rlScalarConstTauAlgInfAndSupport(AC_tau_hl, AV_hl_inf, Foam::GAUR_BATCH_SUPPORT_tau_hl, Foam::GAUR_BATCH_SUPPORT_gInf_hl);

        t[ICaL_d]    = Foam::rlScalarAlgAndSupport(AV_tau_d,    AV_d_inf,    Foam::GAUR_BATCH_SUPPORT_tau_d,    Foam::GAUR_BATCH_SUPPORT_gInf_d);
        t[ICaL_fca]  = Foam::rlScalarAlgAndSupport(AV_tau_fca,  AV_fca_inf,  Foam::GAUR_BATCH_SUPPORT_tau_fca,  Foam::GAUR_BATCH_SUPPORT_gInf_fca);
        t[ICaL_ff]   = Foam::rlScalarAlgAndSupport(AV_tau_ff,   AV_ff_inf,   Foam::GAUR_BATCH_SUPPORT_tau_ff,   Foam::GAUR_BATCH_SUPPORT_gInf_ff);
        t[ICaL_fs]   = Foam::rlScalarAlgAndSupport(AV_tau_fs,   AV_fs_inf,   Foam::GAUR_BATCH_SUPPORT_tau_fs,   Foam::GAUR_BATCH_SUPPORT_gInf_fs);

        t[IKr_xr]    = Foam::rlScalarAlgAndSupport(AV_tau_xr,   AV_xr_inf,   Foam::GAUR_BATCH_SUPPORT_tau_xr,   Foam::GAUR_BATCH_SUPPORT_gInf_xr);

        t[IKs_xs1]   = Foam::rlScalarAlgAndSupport(AV_tau_xs1,  AV_xs1_inf,  Foam::GAUR_BATCH_SUPPORT_tau_xs1,  Foam::GAUR_BATCH_SUPPORT_gInf_xs1);
        t[IKs_xs2]   = Foam::rlScalarAlgAndSupport(AV_tau_xs2,  AV_xs2_inf,  Foam::GAUR_BATCH_SUPPORT_tau_xs2,  Foam::GAUR_BATCH_SUPPORT_gInf_xs2);

        t[CICR_Jrel1] = Foam::rlScalarAlgAndSupport(AV_tau_Jrel1, AV_Jrel1_inf, Foam::GAUR_BATCH_SUPPORT_tau_Jrel1, Foam::GAUR_BATCH_SUPPORT_gInf_Jrel1);
        t[CICR_Jrel2] = Foam::rlScalarAlgAndSupport(AV_tau_Jrel2, AV_Jrel2_inf, Foam::GAUR_BATCH_SUPPORT_tau_Jrel2, Foam::GAUR_BATCH_SUPPORT_gInf_Jrel2);

        return t;
    }();
}



Foam::GaurBatched::GaurBatched
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
        GaurFamilyInfo()
    ),
    useSoAEvaluator_
    (
        dict.lookupOrDefault<Switch>("useSoAEvaluator", false)
    ),
    useCompactSupport_(useGaurCompactSupport(dict)),
    stimulusPOD_()
#ifdef HAS_CUDA
  , useDevice_(false)
#endif
{
    ionicModel::setTissueFromDict();

    if (useSoAEvaluator_ || useCompactSupport_)
    {
        setHotPathSupportSize(NUM_GAUR_BATCH_SUPPORT);
    }

    if (useSoAEvaluator_)
    {
        const word integrator =
            dict.lookupOrDefault<word>("batchedIntegrator", "rushLarsen");
        if (integrator != "euler")
        {
            WarningInFunction
                << "useSoAEvaluator is enabled for " << type()
                << " but `batchedIntegrator " << integrator
                << "` is not honored on the batched path; reverting to "
                << "explicit Euler. Set `batchedIntegrator euler;` to "
                << "silence this warning." << nl;
        }
    }

#ifdef HAS_CUDA
    if (useSoAEvaluator_)
    {
        int nDevices = 0;
        cudaError_t err = cudaGetDeviceCount(&nDevices);
        if (err == cudaSuccess && nDevices > 0)
        {
            const int rank = Pstream::myProcNo();
            CARDIAC_CUDA_CHECK(cudaSetDevice(rank % nDevices));
            useDevice_ = true;
            Info<< "GaurBatched: rank " << rank << " using CUDA device "
                << (rank % nDevices) << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "GaurBatched: useSoAEvaluator is on but no CUDA "
                << "device is visible (" << cudaGetErrorString(err)
                << "); falling back to host SIMD path." << nl;
        }
    }
#endif

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    GaurinitConsts
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


Foam::GaurcompactBatched::GaurcompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    GaurBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

Foam::GaurBatched::~GaurBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}




void Foam::GaurBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_GAUR_BATCH_SUPPORT),
            static_cast<std::size_t>(nCells())
        );
    }
#endif
    configuredBatchedIonicModel::prepareIOAccess(requestedNames, needsAlgebraics);
}


void Foam::GaurBatched::importFields
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


void Foam::GaurBatched::solveODE
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


void Foam::GaurBatched::solveBatched
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
    const bool   solveVm = solveVmWithinODESolver();

    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());

    if (!solveVm)
    {
        for (label cellI = 0; cellI < N; ++cellI)
        {
            state(cellI, cell_v) = vmToState(Vm[cellI]);
        }
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
        const int tFlag = static_cast<int>(tissue());

        cuda_.allocate
        (
            static_cast<std::size_t>(N),
            static_cast<std::size_t>(NUM_STATES),
            static_cast<std::size_t>(NUM_GAUR_BATCH_SUPPORT),
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
            launchGaurBatchKernel
            (
                tSub, cuda_.d_constants, static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                tFlag, solveVm, stimulusPOD_
            );
            launchGaurEulerStepKernel
            (
                cuda_.d_states, cuda_.d_rates,
                static_cast<double>(dtSubstep),
                static_cast<int>(N),
                static_cast<int>(NUM_STATES),
                solveVm,
                static_cast<int>(cell_v)
            );
        }

        launchGaurBatchKernel
        (
            tStart + dtModel, cuda_.d_constants, static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            tFlag, solveVm, stimulusPOD_
        );
        launchGaurScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 1.0,        // Gaur Iion is direct
            static_cast<int>(N),
            static_cast<int>(GAUR_BATCH_SUPPORT_Iion_cm)
        );

        cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));

        cuda_.deviceDirty = true;
        setIOEvaluationModelTime(tStart + dtModel);
        return;
    }
#endif // HAS_CUDA

    for (label sub = 0; sub < nSub; ++sub)
    {
        const scalar tSub = tStart + scalar(sub)*dtSubstep;
        GaurComputeVariablesBatch
        (
            tSub, CONSTS, static_cast<int>(N), 0, static_cast<int>(N),
            STATES_SoA, RATES_SoA, SUPPORT_SoA,
            solveVm, stimulusPOD_
        );

        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            if (!solveVm && stateI == cell_v) continue;
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

    GaurComputeVariablesBatch
    (
        tStart + dtModel, CONSTS,
        static_cast<int>(N), 0, static_cast<int>(N),
        STATES_SoA, RATES_SoA, SUPPORT_SoA,
        solveVm, stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = GAUR_BATCH_SUPPORT_Iion_cm*N;
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = SUPPORT_SoA[IionBase + cellI];
        }
    }

    setIOEvaluationModelTime(tStart + dtModel);
}


Foam::List<Foam::word> Foam::GaurBatched::supportedTissueTypes() const
{
    return {"myocyte"};
}


void Foam::GaurBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    GaurcomputeVariables
    (
        modelTime,
        CONSTANTS_.data(),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        algebraicValues.data(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

Foam::scalar Foam::GaurBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[GAUR_BATCH_SUPPORT_Iion_cm];
}

void Foam::GaurBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (const auto& entry : GaurRushLarsenDispatch)
    {
        projectScalarRushLarsenEntryToSupport
        (
            entry,
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }

    const scalar alpha = algebraics[AV_alpha_aa];
    const scalar beta = algebraics[AV_beta_aa];
    const scalar sum = alpha + beta;
    if (sum > VSMALL && std::isfinite(sum))
    {
        supportValues[GAUR_BATCH_SUPPORT_tau_aa] = 1.0/sum;
        supportValues[GAUR_BATCH_SUPPORT_gInf_aa] = alpha/sum;
    }

    supportValues[GAUR_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}

bool Foam::GaurBatched::rushLarsenParametersFromHotPathSupport
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

    if (stateI == ITo_aa)
    {
        tau = supportValues[GAUR_BATCH_SUPPORT_tau_aa];
        steadyState = supportValues[GAUR_BATCH_SUPPORT_gInf_aa];
        return tau > VSMALL
            && std::isfinite(tau)
            && std::isfinite(steadyState);
    }

    const auto& entry = GaurRushLarsenDispatch[stateI];
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


bool Foam::GaurBatched::rushLarsenParameters
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

    if (stateI == ITo_aa)
    {
        const scalar alpha = algebraicValues[AV_alpha_aa];
        const scalar beta = algebraicValues[AV_beta_aa];
        const scalar sum = alpha + beta;
        if (sum <= VSMALL || !std::isfinite(sum)) return false;

        tau = 1.0/sum;
        steadyState = alpha/sum;
        return std::isfinite(steadyState);
    }

    const auto& entry = GaurRushLarsenDispatch[stateI];
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


void Foam::GaurBatched::derivatives
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
