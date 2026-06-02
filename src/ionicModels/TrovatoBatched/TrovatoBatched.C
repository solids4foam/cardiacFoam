/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "TrovatoBatched.H"
#include "Trovato_2019Batch.H"
#include "batchedRushLarsenEntry.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& TrovatoFamilyInfo();
}

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    TrovatoRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[ICaT_b]       = Foam::rlScalarAlgAndSupport(AV_taub,    AV_bss,    Foam::TROVATO_BATCH_SUPPORT_tau_ICaT_b,    Foam::TROVATO_BATCH_SUPPORT_gInf_ICaT_b);
        t[ICaT_g]       = Foam::rlScalarAlgAndSupport(AV_taug,    AV_gss,    Foam::TROVATO_BATCH_SUPPORT_tau_ICaT_g,    Foam::TROVATO_BATCH_SUPPORT_gInf_ICaT_g);

        t[IK1_xk1]      = Foam::rlScalarAlgAndSupport(AV_txk1,    AV_xk1ss,  Foam::TROVATO_BATCH_SUPPORT_tau_IK1_xk1,   Foam::TROVATO_BATCH_SUPPORT_gInf_IK1_xk1);

        t[IKr_xrf]      = Foam::rlScalarAlgAndSupport(AV_txrf,    AV_xrss,   Foam::TROVATO_BATCH_SUPPORT_tau_IKr_xrf,   Foam::TROVATO_BATCH_SUPPORT_gInf_IKr_xrf);
        t[IKr_xrs]      = Foam::rlScalarAlgAndSupport(AV_txrs,    AV_xrss,   Foam::TROVATO_BATCH_SUPPORT_tau_IKr_xrs,   Foam::TROVATO_BATCH_SUPPORT_gInf_IKr_xrs);

        t[IKs_xs1]      = Foam::rlScalarAlgAndSupport(AV_txs1,    AV_xs1ss,  Foam::TROVATO_BATCH_SUPPORT_tau_IKs_xs1,   Foam::TROVATO_BATCH_SUPPORT_gInf_IKs_xs1);
        t[IKs_xs2]      = Foam::rlScalarAlgAndSupport(AV_txs2,    AV_xs2ss,  Foam::TROVATO_BATCH_SUPPORT_tau_IKs_xs2,   Foam::TROVATO_BATCH_SUPPORT_gInf_IKs_xs2);

        t[INa_hf]       = Foam::rlScalarAlgAndSupport(AV_thf,     AV_hss,    Foam::TROVATO_BATCH_SUPPORT_tau_INa_hf,    Foam::TROVATO_BATCH_SUPPORT_gInf_INa_hf);
        t[INa_hs]       = Foam::rlScalarAlgAndSupport(AV_ths,     AV_hss,    Foam::TROVATO_BATCH_SUPPORT_tau_INa_hs,    Foam::TROVATO_BATCH_SUPPORT_gInf_INa_hs);
        t[INa_hsp]      = Foam::rlScalarAlgAndSupport(AV_thsp,    AV_hssp,   Foam::TROVATO_BATCH_SUPPORT_tau_INa_hsp,   Foam::TROVATO_BATCH_SUPPORT_gInf_INa_hsp);
        t[INa_m]        = Foam::rlScalarAlgAndSupport(AV_tm,      AV_mss,    Foam::TROVATO_BATCH_SUPPORT_tau_INa_m,     Foam::TROVATO_BATCH_SUPPORT_gInf_INa_m);
        t[INa_j]        = Foam::rlScalarAlgAndSupport(AV_tj,      AV_jss,    Foam::TROVATO_BATCH_SUPPORT_tau_INa_j,     Foam::TROVATO_BATCH_SUPPORT_gInf_INa_j);
        t[INa_jp]       = Foam::rlScalarAlgAndSupport(AV_tjp,     AV_jss,    Foam::TROVATO_BATCH_SUPPORT_tau_INa_jp,    Foam::TROVATO_BATCH_SUPPORT_gInf_INa_jp);

        t[If_y]         = Foam::rlScalarAlgAndSupport(AV_tauy,    AV_yss,    Foam::TROVATO_BATCH_SUPPORT_tau_If_y,      Foam::TROVATO_BATCH_SUPPORT_gInf_If_y);

        t[Ito_a]        = Foam::rlScalarAlgAndSupport(AV_taua,    AV_ass,    Foam::TROVATO_BATCH_SUPPORT_tau_Ito_a,     Foam::TROVATO_BATCH_SUPPORT_gInf_Ito_a);
        t[Ito_i1]       = Foam::rlScalarAlgAndSupport(AV_tauis,   AV_iss,    Foam::TROVATO_BATCH_SUPPORT_tau_Ito_i1,    Foam::TROVATO_BATCH_SUPPORT_gInf_Ito_i1);
        t[Ito_i2]       = Foam::rlScalarAlgAndSupport(AV_tauif,   AV_iss,    Foam::TROVATO_BATCH_SUPPORT_tau_Ito_i2,    Foam::TROVATO_BATCH_SUPPORT_gInf_Ito_i2);

        t[INaL_hL]      = Foam::rlScalarConstTauAlgInfAndSupport(AC_thL,   AV_hLss,  Foam::TROVATO_BATCH_SUPPORT_tau_INaL_hL,  Foam::TROVATO_BATCH_SUPPORT_gInf_INaL_hL);
        t[INaL_mL]      = Foam::rlScalarAlgAndSupport(AV_tmL,     AV_mLss,   Foam::TROVATO_BATCH_SUPPORT_tau_INaL_mL,   Foam::TROVATO_BATCH_SUPPORT_gInf_INaL_mL);
        t[INaL_hLp]     = Foam::rlScalarConstTauAlgInfAndSupport(AC_thLp,  AV_hLssp, Foam::TROVATO_BATCH_SUPPORT_tau_INaL_hLp, Foam::TROVATO_BATCH_SUPPORT_gInf_INaL_hLp);

        t[ICaL_d]       = Foam::rlScalarAlgAndSupport(AV_td,      AV_dss,    Foam::TROVATO_BATCH_SUPPORT_tau_ICaL_d,    Foam::TROVATO_BATCH_SUPPORT_gInf_ICaL_d);
        t[ICaL_ff]      = Foam::rlScalarAlgAndSupport(AV_tff,     AV_fss,    Foam::TROVATO_BATCH_SUPPORT_tau_ICaL_ff,   Foam::TROVATO_BATCH_SUPPORT_gInf_ICaL_ff);
        t[ICaL_fs]      = Foam::rlScalarAlgAndSupport(AV_tfs,     AV_fss,    Foam::TROVATO_BATCH_SUPPORT_tau_ICaL_fs,   Foam::TROVATO_BATCH_SUPPORT_gInf_ICaL_fs);
        t[ICaL_fcaf]    = Foam::rlScalarAlgAndSupport(AV_tfcaf,   AV_fcass,  Foam::TROVATO_BATCH_SUPPORT_tau_ICaL_fcaf, Foam::TROVATO_BATCH_SUPPORT_gInf_ICaL_fcaf);
        t[ICaL_fcafp]   = Foam::rlScalarAlgAndSupport(AV_tfcafp,  AV_fcass,  Foam::TROVATO_BATCH_SUPPORT_tau_ICaL_fcafp, Foam::TROVATO_BATCH_SUPPORT_gInf_ICaL_fcafp);
        t[ICaL_fcas]    = Foam::rlScalarAlgAndSupport(AV_tfcas,   AV_fcass,  Foam::TROVATO_BATCH_SUPPORT_tau_ICaL_fcas, Foam::TROVATO_BATCH_SUPPORT_gInf_ICaL_fcas);
        t[ICaL_ffp]     = Foam::rlScalarAlgAndSupport(AV_tffp,    AV_fss,    Foam::TROVATO_BATCH_SUPPORT_tau_ICaL_ffp,  Foam::TROVATO_BATCH_SUPPORT_gInf_ICaL_ffp);
        t[ICaL_jca]     = Foam::rlScalarConstTauAlgInfAndSupport(AC_tjca,  AV_fcass, Foam::TROVATO_BATCH_SUPPORT_tau_ICaL_jca, Foam::TROVATO_BATCH_SUPPORT_gInf_ICaL_jca);

        t[ryr_Jrel1]    = Foam::rlScalarAlgAndSupport(AV_ireltau,  AV_irelss,  Foam::TROVATO_BATCH_SUPPORT_tau_ryr_Jrel1, Foam::TROVATO_BATCH_SUPPORT_gInf_ryr_Jrel1);
        t[ryr_Jrel2]    = Foam::rlScalarAlgAndSupport(AV_ireltau2, AV_irelss2, Foam::TROVATO_BATCH_SUPPORT_tau_ryr_Jrel2, Foam::TROVATO_BATCH_SUPPORT_gInf_ryr_Jrel2);

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
    void launchTrovatoBatchKernel
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

    void launchTrovatoEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchTrovatoRushLarsenStepKernel
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

    void launchTrovatoScaleIonKernel
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
    defineTypeNameAndDebug(TrovatoBatched, 0);
    defineTypeNameAndDebug(TrovatocompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        TrovatocompactBatched,
        dictionary
    );
}

Foam::TrovatoBatched::TrovatoBatched
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
        TrovatoFamilyInfo()
    ),
#ifdef HAS_CUDA
    useDevice_(false),
#endif
    stimulusPOD_()
{
    ionicModel::setTissueFromDict();

    setHotPathSupportSize(NUM_TROVATO_BATCH_SUPPORT);

#ifdef HAS_CUDA
    {
        int nDevices = 0;
        cudaError_t err = cudaGetDeviceCount(&nDevices);
        if (err == cudaSuccess && nDevices > 0)
        {
            const int rank = Pstream::myProcNo();
            const int devId = rank % nDevices;
            CARDIAC_CUDA_CHECK(cudaSetDevice(devId));
            useDevice_ = true;
            Info<< "TrovatoBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
    }
#endif

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    TrovatoinitConsts
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

Foam::TrovatocompactBatched::TrovatocompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    TrovatoBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

Foam::TrovatoBatched::~TrovatoBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}




Foam::List<Foam::word> Foam::TrovatoBatched::supportedTissueTypes() const
{
    return {"myocyte"};
}

void Foam::TrovatoBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    TrovatocomputeVariables
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

Foam::scalar Foam::TrovatoBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[Foam::TROVATO_BATCH_SUPPORT_Iion_cm];
}

void Foam::TrovatoBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (label stateI = 0; stateI < NUM_STATES; ++stateI)
    {
        projectScalarRushLarsenEntryToSupport
        (
            TrovatoRushLarsenDispatch[stateI],
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }
    supportValues[Foam::TROVATO_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}

void Foam::TrovatoBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_TROVATO_BATCH_SUPPORT),
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


void Foam::TrovatoBatched::importFields
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


void Foam::TrovatoBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

void Foam::TrovatocompactBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
#ifdef HAS_CUDA
    if (useDevice_ && !utilitiesMode())
    {
        solveOnDevice(stepStartTime, deltaT, Vm, Im);
        return;
    }
#endif
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

#ifdef HAS_CUDA
void Foam::TrovatocompactBatched::solveOnDevice
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
            state(cellI, membrane_v) = vmToState(Vm[cellI]);
        }
    }

    markIODirty();
    setIOEvaluationModelTime(tStart);

    cuda_.allocate
    (
        static_cast<std::size_t>(N),
        static_cast<std::size_t>(NUM_STATES),
        static_cast<std::size_t>(nHotPathSupport()),
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
        launchTrovatoBatchKernel
        (
            tSub, cuda_.d_constants, static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );
        launchTrovatoRushLarsenStepKernel
        (
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            static_cast<double>(dtSubstep),
            static_cast<int>(N),
            static_cast<int>(NUM_STATES),
            solveVm,
            static_cast<int>(membrane_v)
        );
    }

    launchTrovatoBatchKernel
    (
        tStart + dtModel, cuda_.d_constants, static_cast<int>(N),
        cuda_.d_states, cuda_.d_rates, cuda_.d_support,
        solveVm, stimulusPOD_
    );
    launchTrovatoScaleIonKernel
    (
        cuda_.d_support, cuda_.d_Im, 1.0,
        static_cast<int>(N),
        static_cast<int>(Foam::TROVATO_BATCH_SUPPORT_Iion_cm)
    );
    cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));
    cuda_.deviceDirty = true;
    setIOEvaluationModelTime(tStart + dtModel);
}
#endif // HAS_CUDA


bool Foam::TrovatoBatched::rushLarsenParametersFromHotPathSupport
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

    return resolveSupportRushLarsenEntry
    (
        TrovatoRushLarsenDispatch[stateI],
        CONSTANTS_,
        supportValues,
        VSMALL,
        steadyState,
        tau
    );
}

void Foam::TrovatoBatched::derivatives
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
