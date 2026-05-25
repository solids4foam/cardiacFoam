/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "ToRORd_dynClBatched.H"
#include "ToRORd_dynCl_2023Batch.H"
#include "ToRORd_dynCl_2023Names.H"
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
    defineTypeNameAndDebug(ToRORd_dynClBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        ToRORd_dynClBatched,
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

    if (useSoAEvaluator_)
    {
        setHotPathSupportSize(NUM_TORORD_DYNCL_BATCH_SUPPORT);

        const word integrator =
            dict.lookupOrDefault<word>("batchedIntegrator", "rlGatesHeunRest");
        if (integrator != "euler")
        {
            WarningInFunction
                << "useSoAEvaluator is enabled for "
                << type() << " but `batchedIntegrator " << integrator
                << "` is not yet honored on the batched path; reverting "
                << "to explicit Euler for the SoA solver. Set "
                << "`batchedIntegrator euler;` to silence this warning, "
                << "or unset `useSoAEvaluator` to use the "
                << "validated cell-major Rush-Larsen + Heun path." << nl;
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
            launchToRORd_dynClEulerStepKernel
            (
                cuda_.d_states, cuda_.d_rates,
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




// ************************************************************************* //
