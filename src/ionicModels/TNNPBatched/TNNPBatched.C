/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "TNNPBatched.H"
#include "TNNP_2004Batch.H"
#include "batchedRushLarsenEntry.H"
#include "TNNP_2004Names.H"

namespace Foam
{
    const ionicModelFamilyInfo& TNNPFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            TNNP_CONSTANTS_NAMES,
            TNNP_STATES_NAMES,
            TNNP_ALGEBRAIC_NAMES,
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
#include "clockTime.H"
#include <array>
#include <cstdlib>

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchTnnpBatchKernel
    (
        double t,
        const double* d_CONSTANTS,
        const double* d_CELL_CONSTANTS,
        bool useCellConstants,
        int N,
        const double* d_STATES,
        double* d_RATES,
        double* d_SUPPORT,
        int tissueFlag,
        bool solveVm,
        StimulusProtocolPOD stimulus
    );
    void launchTnnpEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );
    void launchTnnpRushLarsenStepKernel
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
    void launchTnnpScaleIonKernel
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
    defineTypeNameAndDebug(TNNPBatched, 0);
    defineTypeNameAndDebug(TNNPcompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        TNNPcompactBatched,
        dictionary
    );
}

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    TNNPRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[Xr1] = Foam::rlScalarAlgAndSupport(tau_xr1, xr1_inf, Foam::TNNP_BATCH_SUPPORT_tau_xr1, Foam::TNNP_BATCH_SUPPORT_gInf_xr1);
        t[Xr2] = Foam::rlScalarAlgAndSupport(tau_xr2, xr2_inf, Foam::TNNP_BATCH_SUPPORT_tau_xr2, Foam::TNNP_BATCH_SUPPORT_gInf_xr2);
        t[Xs]  = Foam::rlScalarAlgAndSupport(tau_xs,  xs_inf,  Foam::TNNP_BATCH_SUPPORT_tau_xs,  Foam::TNNP_BATCH_SUPPORT_gInf_xs);
        t[m]   = Foam::rlScalarAlgAndSupport(tau_m,   m_inf,   Foam::TNNP_BATCH_SUPPORT_tau_m,   Foam::TNNP_BATCH_SUPPORT_gInf_m);
        t[h]   = Foam::rlScalarAlgAndSupport(tau_h,   h_inf,   Foam::TNNP_BATCH_SUPPORT_tau_h,   Foam::TNNP_BATCH_SUPPORT_gInf_h);
        t[j]   = Foam::rlScalarAlgAndSupport(tau_j,   j_inf,   Foam::TNNP_BATCH_SUPPORT_tau_j,   Foam::TNNP_BATCH_SUPPORT_gInf_j);
        t[d]   = Foam::rlScalarAlgAndSupport(tau_d,   d_inf,   Foam::TNNP_BATCH_SUPPORT_tau_d,   Foam::TNNP_BATCH_SUPPORT_gInf_d);
        t[f]   = Foam::rlScalarAlgAndSupport(tau_f,   f_inf,   Foam::TNNP_BATCH_SUPPORT_tau_f,   Foam::TNNP_BATCH_SUPPORT_gInf_f);
        t[fCa] = Foam::rlNone();   // clamped-rate gate — RL returns false on both paths
        t[s]   = Foam::rlScalarAlgAndSupport(tau_s,   s_inf,   Foam::TNNP_BATCH_SUPPORT_tau_s,   Foam::TNNP_BATCH_SUPPORT_gInf_s);
        t[r]   = Foam::rlScalarAlgAndSupport(tau_r,   r_inf,   Foam::TNNP_BATCH_SUPPORT_tau_r,   Foam::TNNP_BATCH_SUPPORT_gInf_r);
        t[g]   = Foam::rlNone();   // clamped-rate gate — RL returns false on both paths

        return t;
    }();
}

Foam::TNNPBatched::TNNPBatched
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
        TNNPFamilyInfo()
    ),
    stimulusPOD_()
#ifdef HAS_CUDA
  , useDevice_(false)
#endif
{
    ionicModel::setTissueFromDict();
    setHotPathSupportSize(NUM_TNNP_BATCH_SUPPORT);

#ifdef HAS_CUDA
    {
        int nDevices = 0;
        cudaError_t err = cudaGetDeviceCount(&nDevices);
        Info<< "TNNPBatched: cudaGetDeviceCount err="
            << cudaGetErrorString(err)
            << " nDevices=" << nDevices << nl;
        if (err == cudaSuccess && nDevices > 0)
        {
            const int rank = Pstream::myProcNo();
            CARDIAC_CUDA_CHECK(cudaSetDevice(rank % nDevices));
            useDevice_ = true;
            Info<< "TNNPBatched: rank " << rank << " using CUDA device "
                << (rank % nDevices) << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "TNNPBatched: no CUDA device visible ("
                << cudaGetErrorString(err)
                << "); falling back to host path." << nl;
        }
    }
#endif

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    TNNPinitConsts
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

Foam::TNNPBatched::~TNNPBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}


Foam::TNNPcompactBatched::TNNPcompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    TNNPBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}


void Foam::TNNPcompactBatched::solveODE
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
void Foam::TNNPcompactBatched::solveOnDevice
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    const bool reportGpuTimings = std::getenv("CARDIAC_REPORT_GPU_TIMINGS");
    const label N = nCells();
    const label nSub = nSubsteps();
    const scalar dtModel = deltaT*timeScaleFactor();
    const scalar dtSubstep = dtModel/scalar(nSub);
    const scalar tStart = stepStartTime*timeScaleFactor();
    const bool solveVm = solveVmWithinODESolver();
    const int tFlag = static_cast<int>(tissue());
    scalarField vmStateSlice;

    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());

    if (!solveVm)
    {
        vmStateSlice.setSize(N);
        for (label cellI = 0; cellI < N; ++cellI)
        {
            const scalar vmState = vmToState(Vm[cellI]);
            state(cellI, V) = vmState;
            vmStateSlice[cellI] = vmState;
        }
    }

    markIODirty();
    setIOEvaluationModelTime(tStart);

    clockTime stageTimer;

    cuda_.allocate
    (
        static_cast<std::size_t>(N),
        static_cast<std::size_t>(NUM_STATES),
        static_cast<std::size_t>(nHotPathSupport()),
        static_cast<std::size_t>(CONSTANTS_.size())
    );
    const scalar allocateWallTime = stageTimer.timeIncrement();

    cuda_.uploadConstants(CONSTANTS_.cdata(), CONSTANTS_.size());
    scalarField flattenedCellConstants;
    if (hasHeterogeneousConstants())
    {
        flattenHeterogeneousConstants(flattenedCellConstants);
        cuda_.uploadCellConstants
        (
            flattenedCellConstants.cdata(),
            static_cast<std::size_t>(N),
            static_cast<std::size_t>(CONSTANTS_.size())
        );
    }
    const scalar constantsUploadWallTime = stageTimer.timeIncrement();

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
        cuda_.syncStateSliceHostToDevice
        (
            vmStateSlice.cdata(),
            static_cast<std::size_t>(V),
            static_cast<std::size_t>(N)
        );
    }
    const scalar stateUploadWallTime = stageTimer.timeIncrement();

    for (label sub = 0; sub < nSub; ++sub)
    {
        const scalar tSub = tStart + scalar(sub)*dtSubstep;
        launchTnnpBatchKernel
        (
            tSub,
            cuda_.d_constants,
            cuda_.d_cellConstants,
            hasHeterogeneousConstants(),
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            tFlag, solveVm, stimulusPOD_
        );
        launchTnnpRushLarsenStepKernel
        (
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            static_cast<double>(dtSubstep),
            static_cast<int>(N),
            static_cast<int>(NUM_STATES),
            solveVm,
            static_cast<int>(V)
        );
    }
    CARDIAC_CUDA_CHECK(cudaDeviceSynchronize());
    const scalar substepKernelWallTime = stageTimer.timeIncrement();

    launchTnnpBatchKernel
    (
        tStart + dtModel,
        cuda_.d_constants,
        cuda_.d_cellConstants,
        hasHeterogeneousConstants(),
        static_cast<int>(N),
        cuda_.d_states, cuda_.d_rates, cuda_.d_support,
        tFlag, solveVm, stimulusPOD_
    );
    launchTnnpScaleIonKernel
    (
        cuda_.d_support, cuda_.d_Im, 1.0,
        static_cast<int>(N),
        static_cast<int>(TNNP_BATCH_SUPPORT_Iion_cm)
    );
    CARDIAC_CUDA_CHECK(cudaDeviceSynchronize());
    const scalar finalKernelWallTime = stageTimer.timeIncrement();

    cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));
    const scalar imDownloadWallTime = stageTimer.timeIncrement();

    const bool debugGpuIm = std::getenv("CARDIAC_DEBUG_GPU_IM");

    if (debugGpuIm)
    {
        cuda_.syncSupportDeviceToHost
        (
            supportSoAData(),
            static_cast<std::size_t>(nHotPathSupport()),
            static_cast<std::size_t>(nCells())
        );

        const label nReport = (N < 5 ? N : 5);
        scalar minIm = GREAT;
        scalar maxIm = -GREAT;
        scalar sumIm = 0.0;
        scalar minSupportIm = GREAT;
        scalar maxSupportIm = -GREAT;
        scalar sumSupportIm = 0.0;

        for (label cellI = 0; cellI < N; ++cellI)
        {
            const scalar downloadedIm = Im[cellI];
            const scalar supportIm =
                support(cellI, TNNP_BATCH_SUPPORT_Iion_cm);

            minIm = Foam::min(minIm, downloadedIm);
            maxIm = Foam::max(maxIm, downloadedIm);
            sumIm += downloadedIm;

            minSupportIm = Foam::min(minSupportIm, supportIm);
            maxSupportIm = Foam::max(maxSupportIm, supportIm);
            sumSupportIm += supportIm;
        }

        for (label cellI = 0; cellI < nReport; ++cellI)
        {
            Info<< "DEBUG_GPU_IM cell=" << cellI
                << " downloadedIm=" << Im[cellI]
                << " supportIion=" << support(cellI, TNNP_BATCH_SUPPORT_Iion_cm)
                << nl;
        }

        Info<< "DEBUG_GPU_IM stats downloadedIm[min,max,mean]=["
            << minIm << ", " << maxIm << ", " << (sumIm/scalar(N))
            << "] supportIion[min,max,mean]=["
            << minSupportIm << ", " << maxSupportIm << ", "
            << (sumSupportIm/scalar(N)) << "]" << nl;
    }
    const scalar debugSyncWallTime = stageTimer.timeIncrement();

    if (reportGpuTimings)
    {
        Info<< "TNNPBatched GPU timings: cells=" << N
            << " substeps=" << nSub
            << " allocate=" << allocateWallTime << " s"
            << " constantsUpload=" << constantsUploadWallTime << " s"
            << " stateUpload=" << stateUploadWallTime << " s"
            << " substepKernels=" << substepKernelWallTime << " s"
            << " finalKernel=" << finalKernelWallTime << " s"
            << " imDownload=" << imDownloadWallTime << " s";

        if (debugGpuIm)
        {
            Info<< " debugSync=" << debugSyncWallTime << " s";
        }

        Info<< nl;
    }

    cuda_.deviceDirty = true;
    setIOEvaluationModelTime(tStart + dtModel);
}
#endif // HAS_CUDA


void Foam::TNNPBatched::prepareIOAccess
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
            static_cast<std::size_t>(nHotPathSupport()),
            static_cast<std::size_t>(nCells())
        );
    }
#endif
    configuredBatchedIonicModel::prepareIOAccess(requestedNames, needsAlgebraics);
}


void Foam::TNNPBatched::importFields
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


Foam::List<Foam::word> Foam::TNNPBatched::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}


Foam::scalarField Foam::TNNPBatched::constantsForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    TNNPinitConsts
    (
        constants.data(),
        rates.data(),
        states.data(),
        tissueFlag,
        dict()
    );

    ionicModelIO::applyConstantOverrides
    (
        constants,
        TNNP_CONSTANTS_NAMES,
        NUM_CONSTANTS,
        dict(),
        type(),
        tissueFlag
    );

    return constants;
}

void Foam::TNNPBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}



Foam::scalar Foam::TNNPBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[TNNP_BATCH_SUPPORT_Iion_cm];
}

void Foam::TNNPBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    TNNPcomputeRates
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


void Foam::TNNPBatched::evaluateState
(
    const label cellI,
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    scalarField& cellConstants = constants(cellI);

    TNNPcomputeRates
    (
        modelTime,
        cellConstants.data(),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        algebraicValues.data(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

void Foam::TNNPBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (const auto& entry : TNNPRushLarsenDispatch)
    {
        projectScalarRushLarsenEntryToSupport
        (
            entry,
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }
    supportValues[TNNP_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}


void Foam::TNNPBatched::evaluateHotPathState
(
    const label cellI,
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField& cellConstants = constants(cellI);
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(cellI, modelTime, stateValues, rateValues, algebraics);

    for (const auto& entry : TNNPRushLarsenDispatch)
    {
        projectScalarRushLarsenEntryToSupport
        (
            entry,
            cellConstants,
            algebraics,
            supportValues
        );
    }
    supportValues[TNNP_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}


bool Foam::TNNPBatched::rushLarsenParametersFromHotPathSupport
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

    const auto& entry = TNNPRushLarsenDispatch[stateI];
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


bool Foam::TNNPBatched::rushLarsenParametersFromHotPathSupport
(
    const label cellI,
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& supportValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI < 0 || stateI >= NUM_STATES) return false;

    const auto& entry = TNNPRushLarsenDispatch[stateI];
    return resolveSupportRushLarsenEntry
    (
        entry,
        constants(cellI),
        supportValues,
        VSMALL,
        steadyState,
        tau
    );
}

void Foam::TNNPBatched::derivatives
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
