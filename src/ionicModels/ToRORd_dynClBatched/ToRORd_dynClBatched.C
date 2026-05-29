/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "ToRORd_dynClBatched.H"
#include "ToRORd_dynCl_2023Batch.H"
#include "ToRORd_dynCl_2023Names.H"
#include "batchedRushLarsenEntry.H"
#include <array>
namespace Foam
{
    const ionicModelFamilyInfo& ToRORd_dynClFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            ToRORd_dynClCONSTANTS_NAMES,
            ToRORd_dynClSTATES_NAMES,
            ToRORd_dynClALGEBRAIC_NAMES,
            V,
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

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    ToRORd_dynClRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        // INa
        t[INa_m]  = Foam::rlScalarAlgAndSupport(AV_tm,  AV_mss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_INa_m,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_m);
        t[INa_h]  = Foam::rlScalarAlgAndSupport(AV_th,  AV_hss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_INa_h,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_h);
        t[INa_j]  = Foam::rlScalarAlgAndSupport(AV_tj,  AV_jss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_INa_j,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_j);
        t[INa_hp] = Foam::rlScalarAlgAndSupport(AV_th,  AV_hssp, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_INa_hp, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_hp);
        t[INa_jp] = Foam::rlScalarAlgAndSupport(AV_tjp, AV_jss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_INa_jp, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_jp);

        // INaL
        t[INaL_mL]  = Foam::rlScalarAlgAndSupport(AV_tmL, AV_mLss, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_mL,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_mL);
        t[INaL_hL]  = Foam::rlScalarConstTauAlgInfAndSupport(AC_thL,  AV_hLss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_hL,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_hL);
        t[INaL_hLp] = Foam::rlScalarConstTauAlgInfAndSupport(AC_thLp, AV_hLssp, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_hLp, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_hLp);

        // Ito
        t[Ito_a]   = Foam::rlScalarAlgAndSupport(AV_ta,   AV_ass,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_a,   Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_a);
        t[Ito_iF]  = Foam::rlScalarAlgAndSupport(AV_tiF,  AV_iss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iF,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iF);
        t[Ito_iS]  = Foam::rlScalarAlgAndSupport(AV_tiS,  AV_iss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iS,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iS);
        t[Ito_ap]  = Foam::rlScalarAlgAndSupport(AV_ta,   AV_assp, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_ap,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_ap);
        t[Ito_iFp] = Foam::rlScalarAlgAndSupport(AV_tiFp, AV_iss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iFp, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iFp);
        t[Ito_iSp] = Foam::rlScalarAlgAndSupport(AV_tiSp, AV_iss,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iSp, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iSp);

        // ICaL voltage- and Ca-dependent gates (Markov nca_* left as Euler)
        t[ICaL_d]     = Foam::rlScalarAlgAndSupport(AV_td,     AV_dss,   Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_d,     Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_d);
        t[ICaL_ff]    = Foam::rlScalarAlgAndSupport(AV_tff,    AV_fss,   Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_ff,    Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_ff);
        t[ICaL_fs]    = Foam::rlScalarAlgAndSupport(AV_tfs,    AV_fss,   Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fs,    Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fs);
        t[ICaL_fcaf]  = Foam::rlScalarAlgAndSupport(AV_tfcaf,  AV_fcass, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcaf,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcaf);
        t[ICaL_fcas]  = Foam::rlScalarAlgAndSupport(AV_tfcas,  AV_fcass, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcas,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcas);
        t[ICaL_jca]   = Foam::rlScalarConstTauAlgInfAndSupport(AC_tjca,  AV_jcass, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_jca, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_jca);
        t[ICaL_ffp]   = Foam::rlScalarAlgAndSupport(AV_tffp,   AV_fss,   Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_ffp,   Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_ffp);
        t[ICaL_fcafp] = Foam::rlScalarAlgAndSupport(AV_tfcafp, AV_fcass, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcafp, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcafp);

        // IKs
        t[IKs_xs1] = Foam::rlScalarAlgAndSupport(AV_txs1, AV_xs1ss, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_IKs_xs1, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_IKs_xs1);
        t[IKs_xs2] = Foam::rlScalarAlgAndSupport(AV_txs2, AV_xs2ss, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_IKs_xs2, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_IKs_xs2);

        // RyR (HH form, precedent: Trovato ryr_Jrel1/2)
        t[Jrel_np] = Foam::rlScalarAlgAndSupport(AV_tau_rel,  AV_Jrel_inf,  Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_Jrel_np, Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_Jrel_np);
        t[Jrel_p]  = Foam::rlScalarAlgAndSupport(AV_tau_relp, AV_Jrel_infp, Foam::TORORD_DYNCL_BATCH_SUPPORT_tau_Jrel_p,  Foam::TORORD_DYNCL_BATCH_SUPPORT_gInf_Jrel_p);

        // All remaining states (V, CaMKt, concentrations, IKr Markov, ICaL nca)
        // intentionally left at rlNone() => integrated by the executor's
        // explicit-Euler fallback, preserving prior behavior bit-for-bit.

        return t;
    }();
}

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchToRORd_dynClBatchKernel
    (
        double t,
        const double* d_CONSTANTS,
        int N,
        const double* d_STATES,
        double* d_RATES,
        double* d_SUPPORT,
        int tissueFlag,
        bool solveVm,
        StimulusProtocolPOD stimulus
    );

    void launchToRORd_dynClEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchToRORd_dynClRushLarsenStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        const double* d_SUPPORT,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchToRORd_dynClScaleIonKernel
    (
        const double* d_SUPPORT,
        double* d_Im,
        double scale,
        int N,
        int IionBase
    );
}

#endif // HAS_CUDA


namespace Foam
{
    bool useToRORd_dynClCompactSupport(const dictionary& dict)
    {
        const word modelName =
            dict.lookupOrDefault<word>("ionicModel", word::null);

        return modelName == "ToRORd_dynClcompactBatched"
            || dict.lookupOrDefault<Switch>("useCompactSupport", false);
    }

    defineTypeNameAndDebug(ToRORd_dynClBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        ToRORd_dynClBatched,
        dictionary
    );
    defineTypeNameAndDebug(ToRORd_dynClcompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        ToRORd_dynClcompactBatched,
        dictionary
    );
}


Foam::ToRORd_dynClBatched::ToRORd_dynClBatched
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
        ToRORd_dynClFamilyInfo()
    ),
    useSoAEvaluator_
    (
        dict.lookupOrDefault<Switch>("useSoAEvaluator", false)
    ),
    useCompactSupport_(useToRORd_dynClCompactSupport(dict)),
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
            Info<< "ToRORd_dynClBatched: rank " << rank << " using CUDA device "
                << devId << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "ToRORd_dynClBatched: useSoAEvaluator is on but "
                << "no CUDA device is visible (" << cudaGetErrorString(err)
                << "); falling back to the host SIMD path." << nl;
        }
    }
#endif

    if (useSoAEvaluator_ || useCompactSupport_)
    {
        setHotPathSupportSize(NUM_TORORD_DYNCL_BATCH_SUPPORT);
    }

    if (useSoAEvaluator_)
    {
        const word integrator =
            dict.lookupOrDefault<word>("batchedIntegrator", "rushLarsen");
        if (integrator != "euler")
        {
            WarningInFunction
                << "useSoAEvaluator is enabled for "
                << type() << " but `batchedIntegrator " << integrator
                << "` is not yet honored on the batched path; reverting "
                << "to explicit Euler for the SoA solver. Set "
                << "`batchedIntegrator euler;` to silence this warning, "
                << "or unset `useSoAEvaluator` to use the "
                << "validated cell-major Rush-Larsen path." << nl;
        }
    }

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    ToRORd_dynClinitConsts
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


Foam::ToRORd_dynClcompactBatched::ToRORd_dynClcompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ToRORd_dynClBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}


Foam::ToRORd_dynClBatched::~ToRORd_dynClBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}




void Foam::ToRORd_dynClBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_TORORD_DYNCL_BATCH_SUPPORT),
            static_cast<std::size_t>(nCells())
        );
    }
#endif
    configuredBatchedIonicModel::prepareIOAccess(requestedNames, needsAlgebraics);
}


void Foam::ToRORd_dynClBatched::importFields
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


void Foam::ToRORd_dynClBatched::solveODE
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


void Foam::ToRORd_dynClBatched::solveBatched
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
            state(cellI, V) = vmToState(Vm[cellI]);
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
            static_cast<std::size_t>(NUM_TORORD_DYNCL_BATCH_SUPPORT),
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
            launchToRORd_dynClBatchKernel
            (
                tSub, cuda_.d_constants,
                static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                tFlag, solveVm, stimulusPOD_
            );
            launchToRORd_dynClRushLarsenStepKernel
            (
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                static_cast<double>(dtSubstep),
                static_cast<int>(N),
                static_cast<int>(NUM_STATES),
                solveVm,
                static_cast<int>(V)
            );
        }

        launchToRORd_dynClBatchKernel
        (
            tStart + dtModel, cuda_.d_constants,
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            tFlag, solveVm, stimulusPOD_
        );

        launchToRORd_dynClScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 1.0,
            static_cast<int>(N),
            static_cast<int>(TORORD_DYNCL_BATCH_SUPPORT_Iion_cm)
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

        ToRORd_dynClComputeVariablesBatch
        (
            tSub,
            CONSTS,
            static_cast<int>(N),
            0, static_cast<int>(N),
            STATES_SoA,
            RATES_SoA,
            SUPPORT_SoA,
            solveVm,
            stimulusPOD_
        );

        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            if (!solveVm && stateI == V)
            {
                continue;
            }

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
    ToRORd_dynClComputeVariablesBatch
    (
        tEnd,
        CONSTS,
        static_cast<int>(N),
        0, static_cast<int>(N),
        STATES_SoA,
        RATES_SoA,
        SUPPORT_SoA,
        solveVm,
        stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = TORORD_DYNCL_BATCH_SUPPORT_Iion_cm*N;
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


Foam::List<Foam::word> Foam::ToRORd_dynClBatched::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}


void Foam::ToRORd_dynClBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    ToRORd_dynClcomputeVariables
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


Foam::scalar Foam::ToRORd_dynClBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return useCompactSupport_
      ? supportValues[Foam::TORORD_DYNCL_BATCH_SUPPORT_Iion_cm]
      : ionicCurrentFromEvaluation(supportValues);
}


void Foam::ToRORd_dynClBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    if (!useCompactSupport_)
    {
        evaluateState(modelTime, stateValues, rateValues, supportValues);
        return;
    }

    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (label stateI = 0; stateI < NUM_STATES; ++stateI)
    {
        projectScalarRushLarsenEntryToSupport
        (
            ToRORd_dynClRushLarsenDispatch[stateI],
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }
    supportValues[Foam::TORORD_DYNCL_BATCH_SUPPORT_Iion_cm] =
        algebraics[Iion_cm];
}


void Foam::ToRORd_dynClBatched::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(t, y, dydt, algebraics);
}


bool Foam::ToRORd_dynClBatched::rushLarsenParameters
(
    const label stateI,
    const scalarUList& /*stateValues*/,
    const scalarUList& /*rateValues*/,
    const scalarUList& algebraicValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI < 0 || stateI >= NUM_STATES) return false;

    const auto& entry = ToRORd_dynClRushLarsenDispatch[stateI];
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


bool Foam::ToRORd_dynClBatched::rushLarsenParametersFromHotPathSupport
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& supportValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (!useCompactSupport_)
    {
        return rushLarsenParameters
        (
            stateI,
            stateValues,
            rateValues,
            supportValues,
            steadyState,
            tau
        );
    }

    if (stateI < 0 || stateI >= NUM_STATES) return false;

    return resolveSupportRushLarsenEntry
    (
        ToRORd_dynClRushLarsenDispatch[stateI],
        CONSTANTS_,
        supportValues,
        VSMALL,
        steadyState,
        tau
    );
}




// ************************************************************************* //
