/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "TrovatoBatched.H"
#include "Trovato_2019Batch.H"
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& TrovatoFamilyInfo();
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
    addToRunTimeSelectionTable
    (
        ionicModel,
        TrovatoBatched,
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
            Info<< "TrovatoBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "TrovatoBatched: useSoAEvaluator is on but "
                << "no CUDA device is visible (" << cudaGetErrorString(err)
                << "); falling back to the host SIMD path." << nl;
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

    if (useSoAEvaluator_)
    {
        setHotPathSupportSize(NUM_TROVATO_BATCH_SUPPORT);

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
    if (useSoAEvaluator_ && !utilitiesMode())
    {
        solveBatched(stepStartTime, deltaT, Vm, Im);
        return;
    }

    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

void Foam::TrovatoBatched::solveBatched
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
            static_cast<std::size_t>(NUM_TROVATO_BATCH_SUPPORT),
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
                tSub, cuda_.d_constants,
                static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                solveVm, stimulusPOD_
            );

            launchTrovatoEulerStepKernel
            (
                cuda_.d_states, cuda_.d_rates,
                static_cast<double>(dtSubstep),
                static_cast<int>(N),
                static_cast<int>(NUM_STATES),
                solveVm,
                static_cast<int>(membrane_v)
            );
        }

        const scalar tEnd = tStart + dtModel;
        launchTrovatoBatchKernel
        (
            tEnd, cuda_.d_constants,
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );

        launchTrovatoScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 1.0,
            static_cast<int>(N),
            static_cast<int>(TROVATO_BATCH_SUPPORT_Iion_cm)
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

        TrovatoComputeVariablesBatch
        (
            tSub, CONSTS, static_cast<int>(N), 0, static_cast<int>(N),
            STATES_SoA, RATES_SoA, SUPPORT_SoA,
            solveVm, stimulusPOD_
        );

        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            if (!solveVm && stateI == membrane_v) continue;
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
    TrovatoComputeVariablesBatch
    (
        tEnd, CONSTS, static_cast<int>(N), 0, static_cast<int>(N),
        STATES_SoA, RATES_SoA, SUPPORT_SoA,
        solveVm, stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = TROVATO_BATCH_SUPPORT_Iion_cm*N;
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = SUPPORT_SoA[IionBase + cellI];
        }
    }

    setIOEvaluationModelTime(tEnd);
}

bool Foam::TrovatoBatched::rushLarsenParameters
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

    switch (stateI)
    {
        case ICaT_b:
            tau = algebraicValues[AV_taub];
            break;
        case ICaT_g:
            tau = algebraicValues[AV_taug];
            break;
        case IK1_xk1:
            tau = algebraicValues[AV_txk1];
            break;
        case IKr_xrf:
            tau = algebraicValues[AV_txrf];
            break;
        case IKr_xrs:
            tau = algebraicValues[AV_txrs];
            break;
        case IKs_xs1:
            tau = algebraicValues[AV_txs1];
            break;
        case IKs_xs2:
            tau = algebraicValues[AV_txs2];
            break;
        case INa_hf:
            tau = algebraicValues[AV_thf];
            break;
        case INa_hs:
            tau = algebraicValues[AV_ths];
            break;
        case INa_hsp:
            tau = algebraicValues[AV_thsp];
            break;
        case INa_m:
            tau = algebraicValues[AV_tm];
            break;
        case INa_j:
            tau = algebraicValues[AV_tj];
            break;
        case INa_jp:
            tau = algebraicValues[AV_tjp];
            break;
        case If_y:
            tau = algebraicValues[AV_tauy];
            break;
        case Ito_a:
            tau = algebraicValues[AV_taua];
            break;
        case Ito_i1:
            tau = algebraicValues[AV_tauis];
            break;
        case Ito_i2:
            tau = algebraicValues[AV_tauif];
            break;
        case INaL_hL:
            tau = CONSTANTS_[AC_thL];
            break;
        case INaL_mL:
            tau = algebraicValues[AV_tmL];
            break;
        case INaL_hLp:
            tau = CONSTANTS_[AC_thLp];
            break;
        case ICaL_d:
            tau = algebraicValues[AV_td];
            break;
        case ICaL_ff:
            tau = algebraicValues[AV_tff];
            break;
        case ICaL_fs:
            tau = algebraicValues[AV_tfs];
            break;
        case ICaL_fcaf:
            tau = algebraicValues[AV_tfcaf];
            break;
        case ICaL_fcafp:
            tau = algebraicValues[AV_tfcafp];
            break;
        case ICaL_fcas:
            tau = algebraicValues[AV_tfcas];
            break;
        case ICaL_ffp:
            tau = algebraicValues[AV_tffp];
            break;
        case ICaL_jca:
            tau = CONSTANTS_[AC_tjca];
            break;
        case ryr_Jrel1:
            tau = algebraicValues[AV_ireltau];
            break;
        case ryr_Jrel2:
            tau = algebraicValues[AV_ireltau2];
            break;

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
