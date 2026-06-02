/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "CourtemancheBatched.H"
#include "batchedRushLarsenEntry.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& CourtemancheFamilyInfo();
}

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    CourtemancheRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[ina_m]    = Foam::rlScalarAlgAndSupport(AV_ina_m_tau,  AV_ina_m_inf,  Foam::COURTEMANCHE_BATCH_SUPPORT_tau_m,       Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_m);
        t[ina_h]    = Foam::rlScalarAlgAndSupport(AV_ina_h_tau,  AV_ina_h_inf,  Foam::COURTEMANCHE_BATCH_SUPPORT_tau_h,       Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_h);
        t[ina_j]    = Foam::rlScalarAlgAndSupport(AV_ina_j_tau,  AV_ina_j_inf,  Foam::COURTEMANCHE_BATCH_SUPPORT_tau_j,       Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_j);
        t[ical_d]   = Foam::rlScalarAlgAndSupport(AV_ical_d_tau, AV_ical_d_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_d,       Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_d);
        t[ical_f]   = Foam::rlScalarAlgAndSupport(AV_ical_f_tau, AV_ical_f_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_f,       Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_f);
        t[ical_fCa] = Foam::rlScalarConstTauAlgInfAndSupport(AC_ical_fCa_tau, AV_ical_fCa_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_fCa, Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_fCa);
        t[ito_oa]   = Foam::rlScalarAlgAndSupport(AV_ito_oa_tau, AV_ito_oa_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_oa,      Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_oa);
        t[ito_oi]   = Foam::rlScalarAlgAndSupport(AV_ito_oi_tau, AV_ito_oi_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_oi,      Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_oi);
        t[ikur_ua]  = Foam::rlScalarAlgAndSupport(AV_ikur_ua_tau, AV_ikur_ua_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_ua,    Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_ua);
        t[ikur_ui]  = Foam::rlScalarAlgAndSupport(AV_ikur_ui_tau, AV_ikur_ui_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_ui,    Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_ui);
        t[ikr_xr]   = Foam::rlScalarAlgAndSupport(AV_ikr_xr_tau, AV_ikr_xr_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_xr,      Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_xr);
        t[iks_xs]   = Foam::rlScalarAlgAndSupport(AV_iks_xs_tau, AV_iks_xs_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_xs,      Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_xs);
        t[cajsr_u]  = Foam::rlScalarConstTauAlgInfAndSupport(AC_cajsr_u_tau, AV_cajsr_u_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_u, Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_u);
        t[cajsr_w]  = Foam::rlScalarAlgAndSupport(AV_cajsr_w_tau, AV_cajsr_w_inf, Foam::COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_w, Foam::COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_w);

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

    void launchCourtemancheRushLarsenStepKernel
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
    defineTypeNameAndDebug(CourtemanchecompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        CourtemanchecompactBatched,
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
#ifdef HAS_CUDA
    useDevice_(false),
#endif
    stimulusPOD_()
{
    ionicModel::setTissueFromDict();

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
            Info<< "CourtemancheBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
    }
#endif

    setHotPathSupportSize(NUM_COURTEMANCHE_BATCH_SUPPORT);

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

Foam::CourtemanchecompactBatched::CourtemanchecompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    CourtemancheBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

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

Foam::scalar Foam::CourtemancheBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[COURTEMANCHE_BATCH_SUPPORT_Iion_cm];
}

void Foam::CourtemancheBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (const auto& entry : CourtemancheRushLarsenDispatch)
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
        CONSTANTS_[AC_K_rel]
       *stateValues[cajsr_u]*stateValues[cajsr_u]
       *stateValues[cajsr_v]
       *stateValues[cajsr_w]
       *(stateValues[calcium_CaRel] - stateValues[calcium_Cai]);

    const scalar Fn =
        1e-12*CONSTANTS_[AC_V_rel]*Irel
      - 5e-13/CONSTANTS_[AC_F]
       *(0.5*algebraics[AV_ICaL] - 0.2*algebraics[AV_INaCa])
       *CONSTANTS_[AC_Cm];

    supportValues[COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_v] =
        1.0 - 1.0
       /(1.0 + std::exp(-(Fn - 0.2*CONSTANTS_[AC_c1])/CONSTANTS_[AC_c2]));

    supportValues[COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_v] =
        1.91 + 2.09
       /(1.0 + std::exp(-(Fn - CONSTANTS_[AC_c1])/CONSTANTS_[AC_c2]));

    supportValues[COURTEMANCHE_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}

bool Foam::CourtemancheBatched::rushLarsenParametersFromHotPathSupport
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

    if (stateI == cajsr_v)
    {
        tau = supportValues[COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_v];
        steadyState = supportValues[COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_v];
        return tau > VSMALL
            && std::isfinite(tau)
            && std::isfinite(steadyState);
    }

    const auto& entry = CourtemancheRushLarsenDispatch[stateI];
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
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

void Foam::CourtemanchecompactBatched::solveODE
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
void Foam::CourtemanchecompactBatched::solveOnDevice
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
        launchCourtemancheBatchKernel
        (
            tSub, cuda_.d_constants, static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );
        launchCourtemancheRushLarsenStepKernel
        (
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            static_cast<double>(dtSubstep),
            static_cast<int>(N),
            static_cast<int>(NUM_STATES),
            solveVm,
            static_cast<int>(membrane_V)
        );
    }

    launchCourtemancheBatchKernel
    (
        tStart + dtModel, cuda_.d_constants, static_cast<int>(N),
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
    setIOEvaluationModelTime(tStart + dtModel);
}
#endif // HAS_CUDA


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
