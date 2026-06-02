/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "PerisYagueBatched.H"
#include "batchedRushLarsenEntry.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& PerisYagueFamilyInfo();
}

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    PerisYagueRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[ina_m]     = Foam::rlScalarAlgAndSupport(AV_ina_m_tau,    AV_ina_m_inf,    Foam::PERISYAGUE_BATCH_SUPPORT_tau_m,     Foam::PERISYAGUE_BATCH_SUPPORT_gInf_m);
        t[ina_h]     = Foam::rlScalarAlgAndSupport(AV_ina_h_tau,    AV_ina_h_inf,    Foam::PERISYAGUE_BATCH_SUPPORT_tau_h,     Foam::PERISYAGUE_BATCH_SUPPORT_gInf_h);
        t[ina_j]     = Foam::rlScalarAlgAndSupport(AV_ina_j_tau,    AV_ina_j_inf,    Foam::PERISYAGUE_BATCH_SUPPORT_tau_j,     Foam::PERISYAGUE_BATCH_SUPPORT_gInf_j);
        t[ikr_xr]    = Foam::rlScalarAlgAndSupport(AV_ikr_xr_tau,   AV_ikr_xr_inf,   Foam::PERISYAGUE_BATCH_SUPPORT_tau_xr,    Foam::PERISYAGUE_BATCH_SUPPORT_gInf_xr);
        t[iks_xs]    = Foam::rlScalarAlgAndSupport(AV_iks_xs_tau,   AV_iks_xs_inf,   Foam::PERISYAGUE_BATCH_SUPPORT_tau_xs,    Foam::PERISYAGUE_BATCH_SUPPORT_gInf_xs);
        t[ikur_ua]   = Foam::rlScalarAlgAndSupport(AV_ikur_ua_tau,  AV_ikur_ua_inf,  Foam::PERISYAGUE_BATCH_SUPPORT_tau_ua,    Foam::PERISYAGUE_BATCH_SUPPORT_gInf_ua);
        t[ikur_uif]  = Foam::rlScalarAlgAndSupport(AV_ikur_uif_tau, AV_ikur_uif_inf, Foam::PERISYAGUE_BATCH_SUPPORT_tau_uif,   Foam::PERISYAGUE_BATCH_SUPPORT_gInf_uif);
        t[ikur_uis]  = Foam::rlScalarAlgAndSupport(AV_ikur_uis_tau, AV_ikur_uis_inf, Foam::PERISYAGUE_BATCH_SUPPORT_tau_uis,   Foam::PERISYAGUE_BATCH_SUPPORT_gInf_uis);
        t[ical_d]    = Foam::rlScalarAlgAndSupport(AV_ical_d_tau,   AV_ical_d_inf,   Foam::PERISYAGUE_BATCH_SUPPORT_tau_d,     Foam::PERISYAGUE_BATCH_SUPPORT_gInf_d);
        t[ical_f]    = Foam::rlScalarAlgAndSupport(AV_ical_f_tau,   AV_ical_f_inf,   Foam::PERISYAGUE_BATCH_SUPPORT_tau_f,     Foam::PERISYAGUE_BATCH_SUPPORT_gInf_f);
        t[ical_fCa]  = Foam::rlScalarConstTauAlgInfAndSupport(AC_ical_fCa_tau, AV_ical_fCa_inf, Foam::PERISYAGUE_BATCH_SUPPORT_tau_fCa, Foam::PERISYAGUE_BATCH_SUPPORT_gInf_fCa);
        t[iclca_qCa] = Foam::rlScalarLitTauAlgInfAndSupport(2.0, AV_iclca_qCa_inf, Foam::PERISYAGUE_BATCH_SUPPORT_tau_qCa, Foam::PERISYAGUE_BATCH_SUPPORT_gInf_qCa);
        t[ryr_u]     = Foam::rlScalarConstTauAlgInfAndSupport(AC_cajsr_u_tau, AV_ryr_u_inf, Foam::PERISYAGUE_BATCH_SUPPORT_tau_ryr_u, Foam::PERISYAGUE_BATCH_SUPPORT_gInf_ryr_u);
        t[ryr_w]     = Foam::rlScalarAlgAndSupport(AV_ryr_w_tau, AV_ryr_w_inf, Foam::PERISYAGUE_BATCH_SUPPORT_tau_ryr_w, Foam::PERISYAGUE_BATCH_SUPPORT_gInf_ryr_w);

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

    void launchPerisYagueRushLarsenStepKernel
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
    defineTypeNameAndDebug(PerisYaguecompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        PerisYaguecompactBatched,
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
#ifdef HAS_CUDA
    useDevice_(false),
#endif
    stimulusPOD_()
{
    ionicModel::setTissueFromDict();

    setHotPathSupportSize(NUM_PERISYAGUE_BATCH_SUPPORT);

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
            Info<< "PerisYagueBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
    }
#endif

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

    applyIonicConstantOverrides();

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

Foam::PerisYaguecompactBatched::PerisYaguecompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    PerisYagueBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

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

Foam::scalar Foam::PerisYagueBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[PERISYAGUE_BATCH_SUPPORT_Iion_cm];
}

void Foam::PerisYagueBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (const auto& entry : PerisYagueRushLarsenDispatch)
    {
        projectScalarRushLarsenEntryToSupport
        (
            entry,
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }

    const scalar Irel =
        CONSTANTS_[AC_krel]
       *stateValues[ryr_u]*stateValues[ryr_u]
       *stateValues[ryr_v]
       *stateValues[ryr_w]
       *(stateValues[calcium_CaRel] - stateValues[calcium_Cai]);

    const scalar Fn =
        1e-12*CONSTANTS_[AC_V_rel]*Irel
      - 5e-13/CONSTANTS_[AC_F]
       *(0.5*algebraics[AV_ICaL] - 0.2*algebraics[AV_INaCa])
       *CONSTANTS_[AC_Cm];

    supportValues[PERISYAGUE_BATCH_SUPPORT_gInf_ryr_v] =
        1.0 - 1.0
       /(1.0 + std::exp(-(Fn - 0.2*CONSTANTS_[AC_c1])/CONSTANTS_[AC_c2]));

    supportValues[PERISYAGUE_BATCH_SUPPORT_tau_ryr_v] =
        1.91 + 2.09
       /(1.0 + std::exp(-(Fn - CONSTANTS_[AC_c1])/CONSTANTS_[AC_c2]));

    supportValues[PERISYAGUE_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}

bool Foam::PerisYagueBatched::rushLarsenParametersFromHotPathSupport
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

    if (stateI == ryr_v)
    {
        tau = supportValues[PERISYAGUE_BATCH_SUPPORT_tau_ryr_v];
        steadyState = supportValues[PERISYAGUE_BATCH_SUPPORT_gInf_ryr_v];
        return tau > VSMALL
            && std::isfinite(tau)
            && std::isfinite(steadyState);
    }

    const auto& entry = PerisYagueRushLarsenDispatch[stateI];
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
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

void Foam::PerisYaguecompactBatched::solveODE
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
void Foam::PerisYaguecompactBatched::solveOnDevice
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
        launchPerisYagueBatchKernel
        (
            tSub, cuda_.d_constants, static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );
        launchPerisYagueRushLarsenStepKernel
        (
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            static_cast<double>(dtSubstep),
            static_cast<int>(N),
            static_cast<int>(NUM_STATES),
            solveVm,
            static_cast<int>(membrane_V)
        );
    }

    launchPerisYagueBatchKernel
    (
        tStart + dtModel, cuda_.d_constants, static_cast<int>(N),
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
    setIOEvaluationModelTime(tStart + dtModel);
}
#endif // HAS_CUDA


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
