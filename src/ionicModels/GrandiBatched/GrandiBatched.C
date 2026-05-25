/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "GrandiBatched.H"
#include "Grandi_2011Batch.H"
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& GrandiFamilyInfo();
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
    defineTypeNameAndDebug(GrandiBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        GrandiBatched,
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

    if (useSoAEvaluator_)
    {
        setHotPathSupportSize(NUM_GRANDI_BATCH_SUPPORT);

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
    if (stateI < 0 || stateI >= NUM_STATES)
    {
        return false;
    }

    switch (stateI)
    {
        case Ikr_xr:
            tau = algebraicValues[AV_ikr_xr_tau];
            break;
        case Iks_xs:
            tau = algebraicValues[AV_iks_xs_tau];
            break;
        case Ikur_ikur_r:
            tau = algebraicValues[AV_ikur_r_tau];
            break;
        case Ikur_s:
            tau = algebraicValues[AV_ikur_s_tau];
            break;
        case Ina_h:
            tau = algebraicValues[AV_ina_h_tau];
            break;
        case Ina_j:
            tau = algebraicValues[AV_ina_j_tau];
            break;
        case Ina_m:
            tau = algebraicValues[AV_ina_m_tau];
            break;
        case Inal_hl:
            tau = CONSTANTS_[AC_inal_hl_tau];
            break;
        case Inal_ml:
            tau = 1.0/(algebraicValues[AV_inal_ml_a] + algebraicValues[AV_inal_ml_b]);
            break;
        case Ical_d:
            tau = algebraicValues[AV_ical_d_tau];
            break;
        case Ical_f:
            tau = algebraicValues[AV_ical_f_tau];
            break;
        case Ical_fCaB_jn:
            tau = 1.0/(1.7*stateValues[Ca_jn] + 0.0119);
            break;
        case Ical_fCaB_sl:
            tau = 1.0/(1.7*stateValues[Ca_sl] + 0.0119);
            break;
        case Ito_x:
            tau = algebraicValues[AV_ito_x_tau];
            break;
        case Ito_y:
            tau = algebraicValues[AV_ito_y_tau];
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
