/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "AlievPanfilovBatched.H"
#include "AlievPanfilov_1996Batch.H"
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& AlievPanfilovFamilyInfo();
}

#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "stimulusIO.H"

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchAlievPanfilovBatchKernel
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

    void launchAlievPanfilovEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchAlievPanfilovRushLarsenStepKernel
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

    void launchAlievPanfilovScaleIonKernel
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
    defineTypeNameAndDebug(AlievPanfilovBatched, 0);
    defineTypeNameAndDebug(AlievPanfilovcompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        AlievPanfilovcompactBatched,
        dictionary
    );
}

Foam::AlievPanfilovBatched::AlievPanfilovBatched
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
        AlievPanfilovFamilyInfo()
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
            Info<< "AlievPanfilovBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
    }
#endif

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    AlievPanfilovinitConsts
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

    setHotPathSupportSize(NUM_ALIEVPANFILOV_BATCH_SUPPORT);
}

Foam::AlievPanfilovcompactBatched::AlievPanfilovcompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    AlievPanfilovBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

Foam::AlievPanfilovBatched::~AlievPanfilovBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}



Foam::List<Foam::word>
Foam::AlievPanfilovBatched::supportedTissueTypes() const
{
    return {"myocyte"};
}

void Foam::AlievPanfilovBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    AlievPanfilovcomputeVariables
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

Foam::scalar Foam::AlievPanfilovBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return 100.0*supportValues[ALIEVPANFILOV_BATCH_SUPPORT_Iion_cm];
}

void Foam::AlievPanfilovBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    supportValues[ALIEVPANFILOV_BATCH_SUPPORT_tau_recovery_r] =
        1.0/algebraics[AV_eps];
    supportValues[ALIEVPANFILOV_BATCH_SUPPORT_gInf_recovery_r] =
        -CONSTANTS_[AC_k]*stateValues[u]
       *(stateValues[u] - (CONSTANTS_[AC_a] + 1.0));
    supportValues[ALIEVPANFILOV_BATCH_SUPPORT_Iion_cm] =
        algebraics[Iion_cm];
}

void Foam::AlievPanfilovBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_ALIEVPANFILOV_BATCH_SUPPORT),
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


void Foam::AlievPanfilovBatched::importFields
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


void Foam::AlievPanfilovBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

void Foam::AlievPanfilovcompactBatched::solveODE
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
void Foam::AlievPanfilovcompactBatched::solveOnDevice
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
            state(cellI, u) = vmToState(Vm[cellI]);
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
        launchAlievPanfilovBatchKernel
        (
            tSub, cuda_.d_constants, static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );
        launchAlievPanfilovRushLarsenStepKernel
        (
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            static_cast<double>(dtSubstep),
            static_cast<int>(N),
            static_cast<int>(NUM_STATES),
            solveVm,
            static_cast<int>(u)
        );
    }

    launchAlievPanfilovBatchKernel
    (
        tStart + dtModel, cuda_.d_constants, static_cast<int>(N),
        cuda_.d_states, cuda_.d_rates, cuda_.d_support,
        solveVm, stimulusPOD_
    );
    launchAlievPanfilovScaleIonKernel
    (
        cuda_.d_support, cuda_.d_Im, 100.0,
        static_cast<int>(N),
        static_cast<int>(ALIEVPANFILOV_BATCH_SUPPORT_Iion_cm)
    );
    cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));
    cuda_.deviceDirty = true;
    setIOEvaluationModelTime(tStart + dtModel);
}
#endif // HAS_CUDA


bool Foam::AlievPanfilovBatched::rushLarsenParametersFromHotPathSupport
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& supportValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI != recovery_r) return false;

    tau = supportValues[ALIEVPANFILOV_BATCH_SUPPORT_tau_recovery_r];
    if (tau <= VSMALL || !std::isfinite(tau)) return false;

    steadyState = supportValues[ALIEVPANFILOV_BATCH_SUPPORT_gInf_recovery_r];
    return std::isfinite(steadyState);
}

void Foam::AlievPanfilovBatched::derivatives
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
