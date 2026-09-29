/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "monodomainFDAManufacturedBatched.H"
#include "ionicModelIO.H"
#include "ionicSelector.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchMonodomainFDAManufacturedBatchKernel
    (
        const double* d_CONSTANTS,
        int N,
        const double* d_STATES,
        double* d_RATES,
        double* d_SUPPORT
    );
    void launchMonodomainFDAManufacturedEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );
    void launchMonodomainFDAManufacturedScaleIonKernel
    (
        const double* d_SUPPORT,
        double* d_Im,
        double scale,
        int N,
        int IionSlot
    );
    void launchMonodomainFDAManufacturedSetVmKernel
    (
        double* d_STATES,
        const double* d_VM_START,
        const double* d_VM_RATE,
        double modelTimeOffset,
        int N,
        int vmStateI
    );
}
#endif

namespace Foam
{
namespace
{
    const ionicModelFamilyInfo& monodomainFDAManufacturedBatchInfo()
    {
        static const char* const constantNames[NUM_CONSTANTS] =
        {
            "Cm", "Beta", "Chi"
        };
        static const char* const stateNames[NUM_STATES] =
        {
            "V", "u1", "u2", "u3"
        };
        static const char* const algebraicNames[NUM_ALGEBRAIC] =
        {
            "Iion"
        };
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            constantNames,
            stateNames,
            algebraicNames,
            V,
            1.0,
            1.0,
            0.0,
            nullptr
        };
        return info;
    }
}

defineTypeNameAndDebug(monodomainFDAManufacturedBatched, 0);
addToRunTimeSelectionTable
(
    ionicModel,
    monodomainFDAManufacturedBatched,
    dictionary
);


Foam::monodomainFDAManufacturedBatched::monodomainFDAManufacturedBatched
(
    const dictionary& dict,
    const label nIntegrationPoints,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    configuredBatchedIonicModel
    (
        dict,
        nIntegrationPoints,
        initialDeltaT,
        solveVmWithinODESolver,
        NUM_STATES,
        NUM_ALGEBRAIC,
        monodomainFDAManufacturedBatchInfo()
    ),
    manufacturedSourceTermPtr_(nullptr),
    manufacturedSourceChi_(1.0),
    manufacturedSourceCm_(1.0),
    manufacturedSourceStart_(0)
#ifdef HAS_CUDA
  , useDevice_(false)
  , dVmStateStart_(nullptr)
  , dVmStateRate_(nullptr)
  , vmBufferSize_(0)
  , vmStateStartHost_()
  , vmStateRateHost_()
#endif
{
    setTissue(ionicSelector::selectDimension(dict, supportedDimensions()));
    CONSTANTS_ = constantsForTissue(tissue());
    applyIonicConstantOverrides();
    setHotPathSupportSize(NUM_MONODOMAIN_FDA_BATCH_SUPPORT);

    const scalarField initialStates = initialStatesForTissue(tissue());
    for (label cellI = 0; cellI < nCells(); ++cellI)
    {
        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            state(cellI, stateI) = initialStates[stateI];
            rate(cellI, stateI) = 0.0;
        }
    }

    configurePersistentAlgebraics();
    syncAllToIO();

#ifdef HAS_CUDA
    int nDevices = 0;
    const cudaError_t err = cudaGetDeviceCount(&nDevices);
    if (err == cudaSuccess && nDevices > 0)
    {
        const int rank = Pstream::myProcNo();
        const int device = rank % nDevices;
        CARDIAC_CUDA_CHECK(cudaSetDevice(device));
        useDevice_ = true;
        Info<< "monodomainFDAManufacturedBatched: rank " << rank
            << " using CUDA device " << device << " of " << nDevices << nl;
    }
    else
    {
        WarningInFunction
            << "monodomainFDAManufacturedBatched: no CUDA device visible ("
            << cudaGetErrorString(err) << "); falling back to host path."
            << nl;
    }
#endif
}


Foam::monodomainFDAManufacturedBatched::~monodomainFDAManufacturedBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
        if (dVmStateStart_)
        {
            cudaFree(dVmStateStart_);
            dVmStateStart_ = nullptr;
        }
        if (dVmStateRate_)
        {
            cudaFree(dVmStateRate_);
            dVmStateRate_ = nullptr;
        }
        vmBufferSize_ = 0;
    }
#endif
}


Foam::scalarField
Foam::monodomainFDAManufacturedBatched::constantsForTissue
(
    const label dimension
) const
{
    if (dimension < 1 || dimension > 3)
    {
        FatalErrorInFunction
            << "Invalid manufactured monodomain dimension " << dimension
            << exit(FatalError);
    }

    scalarField values(NUM_CONSTANTS, 0.0);
    values[Cm] = 2.0;
    values[Chi] = 3.0;
    values[Beta] = dimension == 1 ? -1.1 : dimension == 2 ? -5.9 : -8.6;
    return values;
}


Foam::scalarField
Foam::monodomainFDAManufacturedBatched::initialStatesForTissue
(
    const label dimension
) const
{
    (void)dimension;
    return scalarField(NUM_STATES, 0.0);
}


void Foam::monodomainFDAManufacturedBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    (void)modelTime;
    const scalar v = stateValues[V];
    const scalar s1 = stateValues[u1];
    const scalar s2 = stateValues[u2];
    const scalar s3 = stateValues[u3];
    const scalar q = s1 + s3 - v;
    const scalar s2sq = s2*s2;
    const scalar vMinusU3 = v - s3;

    rateValues[V] = 0.0;
    rateValues[u1] = q*q*s2sq + 0.5*q*s2sq*vMinusU3;
    rateValues[u2] = -q*s2sq*s2;
    rateValues[u3] = 0.0;
    algebraicValues[Iion] =
        -0.5*CONSTANTS_[Cm]*q*s2sq*vMinusU3
      + (CONSTANTS_[Beta]*vMinusU3)/CONSTANTS_[Chi];
}


void Foam::monodomainFDAManufacturedBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);
    supportValues[MONODOMAIN_FDA_BATCH_SUPPORT_Iion] = algebraics[Iion];
}


Foam::scalar
Foam::monodomainFDAManufacturedBatched::manufacturedSourceCorrection
(
    const label cellI
) const
{
    if (!manufacturedSourceTermPtr_)
    {
        return 0.0;
    }

    const label sourceI = manufacturedSourceStart_ + cellI;
    if (sourceI < 0 || sourceI >= manufacturedSourceTermPtr_->size())
    {
        FatalErrorInFunction
            << "Manufactured source index " << sourceI
            << " is outside source field size "
            << manufacturedSourceTermPtr_->size()
            << exit(FatalError);
    }

    return (*manufacturedSourceTermPtr_)[sourceI]
         / (manufacturedSourceChi_*manufacturedSourceCm_);
}


void Foam::monodomainFDAManufacturedBatched::setManufacturedSourceTerm
(
    const scalarField& sourceTerm,
    scalar chi,
    scalar CmValue,
    label sourceStart
)
{
    if (mag(chi*CmValue) <= VSMALL)
    {
        FatalErrorInFunction
            << "Manufactured source correction requires non-zero chi*Cm."
            << exit(FatalError);
    }

    manufacturedSourceTermPtr_ = &sourceTerm;
    manufacturedSourceChi_ = chi;
    manufacturedSourceCm_ = CmValue;
    manufacturedSourceStart_ = sourceStart;
}


void Foam::monodomainFDAManufacturedBatched::clearManufacturedSourceTerm()
{
    manufacturedSourceTermPtr_ = nullptr;
    manufacturedSourceChi_ = 1.0;
    manufacturedSourceCm_ = 1.0;
    manufacturedSourceStart_ = 0;
}


void Foam::monodomainFDAManufacturedBatched::solveODE
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
        for (label cellI = 0; cellI < nCells(); ++cellI)
        {
            Im[cellI] += manufacturedSourceCorrection(cellI);
        }
        return;
    }
#endif

    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
    for (label cellI = 0; cellI < nCells(); ++cellI)
    {
        Im[cellI] += manufacturedSourceCorrection(cellI);
    }
}


#ifdef HAS_CUDA
void Foam::monodomainFDAManufacturedBatched::solveOnDevice
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
    const bool extrapolateVm = !solveVm && hasVmRate();

    if (extrapolateVm)
    {
        vmStateStartHost_.setSize(N);
        vmStateRateHost_.setSize(N);
        for (label cellI = 0; cellI < N; ++cellI)
        {
            const scalar vmStart = vmToState(Vm[cellI]);
            const scalar vmEnd =
                vmToState(Vm[cellI] + VmRateSI(cellI)*deltaT);
            vmStateStartHost_[cellI] = vmStart;
            vmStateRateHost_[cellI] =
                (mag(dtModel) > VSMALL)
              ? (vmEnd - vmStart)/dtModel
              : 0.0;
        }
        if (vmBufferSize_ != static_cast<std::size_t>(N))
        {
            if (dVmStateStart_)
            {
                CARDIAC_CUDA_CHECK(cudaFree(dVmStateStart_));
                dVmStateStart_ = nullptr;
            }
            if (dVmStateRate_)
            {
                CARDIAC_CUDA_CHECK(cudaFree(dVmStateRate_));
                dVmStateRate_ = nullptr;
            }

            CARDIAC_CUDA_CHECK
            (
                cudaMalloc
                (
                    reinterpret_cast<void**>(&dVmStateStart_),
                    sizeof(double)*static_cast<std::size_t>(N)
                )
            );
            CARDIAC_CUDA_CHECK
            (
                cudaMalloc
                (
                    reinterpret_cast<void**>(&dVmStateRate_),
                    sizeof(double)*static_cast<std::size_t>(N)
                )
            );
            vmBufferSize_ = static_cast<std::size_t>(N);
        }

        CARDIAC_CUDA_CHECK
        (
            cudaMemcpy
            (
                dVmStateStart_, vmStateStartHost_.cdata(),
                sizeof(double)*static_cast<std::size_t>(N),
                cudaMemcpyHostToDevice
            )
        );
        CARDIAC_CUDA_CHECK
        (
            cudaMemcpy
            (
                dVmStateRate_, vmStateRateHost_.cdata(),
                sizeof(double)*static_cast<std::size_t>(N),
                cudaMemcpyHostToDevice
            )
        );
    }

    if (!solveVm)
    {
        for (label cellI = 0; cellI < N; ++cellI)
        {
            state(cellI, V) = vmToState(Vm[cellI]);
        }
    }

    markIODirty();
    setIOEvaluationModelTime(tStart);
    cuda_.allocate
    (
        static_cast<std::size_t>(N),
        static_cast<std::size_t>(NUM_STATES),
        static_cast<std::size_t>(NUM_MONODOMAIN_FDA_BATCH_SUPPORT),
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
    else if (!solveVm)
    {
        const std::size_t vmStateI = static_cast<std::size_t>(V);
        cuda_.syncStateSliceHostToDevice
        (
            statesSoAData() + vmStateI*static_cast<std::size_t>(N),
            vmStateI,
            static_cast<std::size_t>(N)
        );
    }

    for (label sub = 0; sub < nSub; ++sub)
    {
        if (extrapolateVm)
        {
            setDeviceVmExtrapolant(scalar(sub)*dtSubstep, N);
        }

        launchMonodomainFDAManufacturedBatchKernel
        (
            cuda_.d_constants,
            static_cast<int>(N),
            cuda_.d_states,
            cuda_.d_rates,
            cuda_.d_support
        );
        launchMonodomainFDAManufacturedEulerStepKernel
        (
            cuda_.d_states,
            cuda_.d_rates,
            static_cast<double>(dtSubstep),
            static_cast<int>(N),
            static_cast<int>(NUM_STATES),
            solveVm,
            static_cast<int>(V)
        );
    }

    if (extrapolateVm)
    {
        setDeviceVmExtrapolant(dtModel, N);
    }
    launchMonodomainFDAManufacturedBatchKernel
    (
        cuda_.d_constants,
        static_cast<int>(N),
        cuda_.d_states,
        cuda_.d_rates,
        cuda_.d_support
    );
    launchMonodomainFDAManufacturedScaleIonKernel
    (
        cuda_.d_support,
        cuda_.d_Im,
        1.0/CONSTANTS_[Cm],
        static_cast<int>(N),
        static_cast<int>(MONODOMAIN_FDA_BATCH_SUPPORT_Iion)
    );
    cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));
    cuda_.deviceDirty = true;
    setIOEvaluationModelTime(tStart + dtModel);
}


void Foam::monodomainFDAManufacturedBatched::setDeviceVmExtrapolant
(
    const double modelTimeOffset,
    const label N
)
{
    launchMonodomainFDAManufacturedSetVmKernel
    (
        cuda_.d_states,
        dVmStateStart_,
        dVmStateRate_,
        modelTimeOffset,
        static_cast<int>(N),
        static_cast<int>(V)
    );
}
#endif


void Foam::monodomainFDAManufacturedBatched::evaluateIonicCurrent
(
    const scalar t,
    const scalarField& Vm,
    scalarField& Im
)
{
    configuredBatchedIonicModel::evaluateIonicCurrent(t, Vm, Im);
    for (label cellI = 0; cellI < nCells(); ++cellI)
    {
        Im[cellI] += manufacturedSourceCorrection(cellI);
    }
}


void Foam::monodomainFDAManufacturedBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_MONODOMAIN_FDA_BATCH_SUPPORT),
            static_cast<std::size_t>(nCells())
        );
    }
#endif
    configuredBatchedIonicModel::prepareIOAccess(requestedNames, needsAlgebraics);
}


void Foam::monodomainFDAManufacturedBatched::importFields
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


void Foam::monodomainFDAManufacturedBatched::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(t, y, dydt, algebraics);
    if (!solveVmWithinODESolver())
    {
        dydt[V] = activeVmRate();
    }
}

} // End namespace Foam
