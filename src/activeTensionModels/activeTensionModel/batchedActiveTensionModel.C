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

Class
    Foam::batchedActiveTensionModel

Description
    Reusable base for GPU-friendly active tension models that execute as batched,
    explicit cell updates while preserving the existing activeTensionModel API.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "batchedActiveTensionModel.H"
#include "restartStateIO.H"
#ifdef HAS_CUDA
#include "Pstream.H"
#include <chrono>
#endif

namespace Foam
{

batchedActiveTensionModel::batchedActiveTensionModel
(
    const dictionary& dict,
    const label nIntegrationPoints,
    const label nStates,
    const label nAlgebraics
)
:
    activeTensionModel(dict, nIntegrationPoints),
    nCells_(nIntegrationPoints),
    nStates_(nStates),
    nAlgebraics_(nAlgebraics),
    nSubsteps_(dict.lookupOrDefault<label>("batchedSubsteps", 1)),
    persistAlgebraics_(dict.lookupOrDefault<Switch>("storeBatchedAlgebraics", true)),
    parallelCellUpdates_(dict.lookupOrDefault<Switch>("batchedParallelCells", false)),
    parallelMinCells_(dict.lookupOrDefault<label>("batchedParallelMinCells", 256)),
    core_(nIntegrationPoints, nStates, nAlgebraics),
    ioStates_(nIntegrationPoints),
    ioRates_(nIntegrationPoints),
    ioAlgebraics_(nIntegrationPoints),
    ioSynchronized_(false),
    solveScratch_(),
    solveScratchThreads_(-1)
#ifdef HAS_CUDA
    ,
    cuda_(),
    useCUDA_(dict.lookupOrDefault<Switch>("batchedUseCUDA", true)),
    cudaProfile_(dict.lookupOrDefault<Switch>("batchedCUDAProfile", false)),
    cudaEnabled_(false),
    cudaAnnounced_(false),
    cudaHostStateStale_(false),
    cudaProfileCalls_(0),
    cudaProfileH2DSeconds_(0.0),
    cudaProfileKernelSeconds_(0.0),
    cudaProfileD2HSeconds_(0.0),
    cudaProfileTotalSeconds_(0.0)
#endif
{
    if (nSubsteps_ < 1) nSubsteps_ = 1;
    if (parallelMinCells_ < 1) parallelMinCells_ = 1;

#ifndef _OPENMP
    if (parallelCellUpdates_)
    {
        WarningInFunction
            << "batchedParallelCells was requested for active tension model "
            << type()
            << " but this library was built without OpenMP support. "
            << "Falling back to serial cell loops." << nl;
        parallelCellUpdates_ = false;
    }
#endif

    core_.setAlgebraicsStorageEnabled(persistAlgebraics_);

#ifdef HAS_CUDA
    // CUDA evaluation needs an algebraic mirror even when callers elect not
    // to expose/store algebraics through the restart/output API.
    if (!persistAlgebraics_)
    {
        core_.setAlgebraicsStorageEnabled(true);
    }

    int nDevices = 0;
    const cudaError_t cudaStatus = cudaGetDeviceCount(&nDevices);
    if (useCUDA_ && cudaStatus == cudaSuccess && nDevices > 0)
    {
        const int device = Pstream::myProcNo() % nDevices;
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaSetDevice(device));
        cudaEnabled_ = true;
    }
#endif

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        ioStates_.set(cellI, new scalarField(nStates_, 0.0));
        ioRates_.set(cellI, new scalarField(nStates_, 0.0));

        if (persistAlgebraics_)
        {
            ioAlgebraics_.set(cellI, new scalarField(nAlgebraics_, 0.0));
        }
    }
}


#ifdef HAS_CUDA
bool batchedActiveTensionModel::solveOnCUDA
(
    const scalar t,
    const scalar dt,
    const scalarField& lambda,
    scalarField& Ta
) const
{
    if (!cudaEnabled_)
    {
        return false;
    }

    const auto totalStart = std::chrono::steady_clock::now();
    cudaEvent_t timingStart = nullptr;
    cudaEvent_t timingStop = nullptr;
    if (cudaProfile_)
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventCreate(&timingStart));
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventCreate(&timingStop));
    }
    const auto accumulateElapsed =
        [&](double& accumulator)
        {
            if (!cudaProfile_) return;
            CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventRecord(timingStop));
            CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventSynchronize(timingStop));
            float milliseconds = 0.0f;
            CARDIAC_ACTIVE_TENSION_CUDA_CHECK
            (
                cudaEventElapsedTime(&milliseconds, timingStart, timingStop)
            );
            accumulator += 1.0e-3*milliseconds;
        };

    const scalarField* constants = ioConstantsPtr();
    if (!constants)
    {
        FatalErrorInFunction
            << "CUDA active-tension model " << type()
            << " did not supply its constants." << abort(FatalError);
    }

    cuda_.allocate(nCells_, nStates_, nAlgebraics_, constants->size());
    if (!cuda_.constantsUploaded)
    {
        cuda_.uploadConstants(constants->cdata(), constants->size());
    }
    if (!cuda_.stateResident)
    {
        cuda_.uploadInitialFields
        (
            core_.statesSoAData(), core_.algebraicsSoAData(),
            driveSignals_.data(), nCells_, nStates_, nAlgebraics_
        );
        cuda_.stateResident = true;
    }
    prepareCUDAInputs(driveSignals_, lambda);
    if (cudaProfile_)
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventRecord(timingStart));
    }
    uploadCUDAInputs();
    accumulateElapsed(cudaProfileH2DSeconds_);

    const scalar dtModel = dt*timeScaleFactor()/scalar(nSubsteps_);
    if (cudaProfile_)
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventRecord(timingStart));
    }
    for (label substep = 0; substep < nSubsteps_; ++substep)
    {
        launchCUDAKernel();
        launchBatchedActiveTensionEulerKernel
        (
            cuda_.d_states, cuda_.d_rates, nCells_, nStates_, dtModel
        );
    }

    // Preserve the host executor's post-step convention: output rates and
    // algebraics are evaluated at the advanced state, not at the last Euler
    // right-hand-side state.
    launchCUDAKernel();
    accumulateElapsed(cudaProfileKernelSeconds_);
    if (cudaProfile_)
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventRecord(timingStart));
    }
    downloadCUDATension(Ta);
    accumulateElapsed(cudaProfileD2HSeconds_);
    cudaHostStateStale_ = true;

    if (cudaProfile_)
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventDestroy(timingStart));
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaEventDestroy(timingStop));
        ++cudaProfileCalls_;
        cudaProfileTotalSeconds_ += std::chrono::duration<double>
        (
            std::chrono::steady_clock::now() - totalStart
        ).count();
    }

    if (!cudaAnnounced_)
    {
        Info<< type() << ": rank " << Pstream::myProcNo()
            << " using CUDA device for active-tension ODE updates" << nl;
        cudaAnnounced_ = true;
    }

    (void)t;
    return true;
}
#endif


batchedActiveTensionModel::~batchedActiveTensionModel()
{
#ifdef HAS_CUDA
    if (cudaProfile_ && cudaProfileCalls_ > 0)
    {
        Info<< type() << " CUDA profile: calls=" << cudaProfileCalls_
            << " H2D_s=" << cudaProfileH2DSeconds_
            << " kernels_s=" << cudaProfileKernelSeconds_
            << " D2H_s=" << cudaProfileD2HSeconds_
            << " total_s=" << cudaProfileTotalSeconds_ << nl;
    }
#endif
}


void batchedActiveTensionModel::prepareSolveScratch(const label nThreads) const
{
    if (solveScratchThreads_ == nThreads && solveScratch_.size() == nThreads)
    {
        return;
    }

    solveScratch_.clear();
    solveScratch_.setSize(nThreads);

    for (label threadI = 0; threadI < nThreads; ++threadI)
    {
        solveScratch_.set
        (
            threadI,
            new CellScratch(nStates_, nAlgebraics_)
        );
    }

    solveScratchThreads_ = nThreads;
}


void batchedActiveTensionModel::syncAllToIO() const
{
    if (ioSynchronized_) return;

#ifdef HAS_CUDA
    if (cudaHostStateStale_)
    {
        cuda_.downloadOutputs
        (
            core_.statesSoAData(), core_.ratesSoAData(), core_.algebraicsSoAData(),
            nCells_, nStates_, nAlgebraics_
        );
        cudaHostStateStale_ = false;
    }
#endif

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        scalarField& stateIO = ioStates_[cellI];
        scalarField& rateIO = ioRates_[cellI];

        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            stateIO[stateI] = state(cellI, stateI);
            rateIO[stateI] = rate(cellI, stateI);
        }

        if (persistAlgebraics_)
        {
            scalarField& algIO = ioAlgebraics_[cellI];
            for (label algI = 0; algI < nAlgebraics_; ++algI)
            {
                algIO[algI] = algebraic(cellI, algI);
            }
        }
    }

    ioSynchronized_ = true;
}


void batchedActiveTensionModel::syncStatesFromIO() const
{
    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        const scalarField& stateIO = ioStates_[cellI];
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            state(cellI, stateI) = stateIO[stateI];
        }
    }
    ioSynchronized_ = false;
#ifdef HAS_CUDA
    // A restart/pre-pacing/state injection changed the host source of truth.
    cuda_.stateResident = false;
    cudaHostStateStale_ = false;
#endif
}


bool batchedActiveTensionModel::supportsRestartState() const
{
    return true;
}


bool batchedActiveTensionModel::readRestartState(const fvMesh& mesh)
{
    syncAllToIO();
    if (!restartStateIO::readStates(mesh, type(), nStates_, ioStates_))
    {
        return false;
    }

    syncStatesFromIO();
    core_.clearTransientSolveData(persistAlgebraics_);
    return true;
}


bool batchedActiveTensionModel::setStates(const UList<scalarField>& states)
{
    syncAllToIO();

    if (ioStates_.size() != states.size())
    {
        return false;
    }

    forAll(states, i)
    {
        ioStates_[i] = states[i];
    }

    syncStatesFromIO();
    core_.clearTransientSolveData(persistAlgebraics_);
    return true;
}


void batchedActiveTensionModel::writeRestartState(const fvMesh& mesh) const
{
    syncAllToIO();
    restartStateIO::writeStates(mesh, type(), nStates_, ioStates_);
}


void batchedActiveTensionModel::refreshRestartState(const fvMesh& mesh)
{
    syncAllToIO();
    CellScratch scratch(nStates_, nAlgebraics_);
    BatchedTensionBackend backend(*this);

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        backend.gatherCellState(cellI, scratch.stateValues);
        scratch.resetPrimary();
        backend.evaluateScratchAtTime
        (
            cellI,
            mesh.time().value()*timeScaleFactor(),
            coupledDriveSignal(cellI),
            1.0,
            scratch
        );
        backend.syncEvaluatedOutputs(cellI, scratch);
    }

    ioSynchronized_ = false;
}


bool batchedActiveTensionModel::restartTension(scalarField& Ta) const
{
    if (Ta.size() != nCells_)
    {
        return false;
    }

    CellScratch scratch(nStates_, nAlgebraics_);
    BatchedTensionBackend backend(*this);

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        backend.gatherCellState(cellI, scratch.stateValues);
        scratch.resetPrimary();
        backend.evaluateScratchAtTime
        (
            cellI,
            0.0,
            coupledDriveSignal(cellI),
            1.0,
            scratch
        );
        Ta[cellI] = backend.activeTensionFromScratch(cellI, scratch);
    }

    return true;
}


void batchedActiveTensionModel::calculateTension
(
    const scalar t,
    const scalar dt,
    const scalarField& lambda,
    scalarField& Ta
)
{
    if (Ta.size() != nCells_)
    {
        FatalErrorInFunction
            << "Ta.size() (" << Ta.size()
            << ") != nCells (" << nCells_ << ")"
            << abort(FatalError);
    }

    currentT_  = t;
    currentDt_ = dt;

    // First time synchronisation from original non-batched IO if needed
    // Assuming initial states were set into core_ or we sync from ioStates_
    // For safety, we trust core_ is the source of truth if we use batched.

    if (driveSignals_.size() != nCells_)
    {
        driveSignals_.setSize(nCells_);
    }

    // Gather per-cell drive signals safely before parallel/GPU dispatch
    // This avoids calling virtual functions inside OpenMP/CUDA kernels
    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        driveSignals_[cellI] = coupledDriveSignal(cellI);
    }

#ifdef HAS_CUDA
    if (solveOnCUDA(t, dt, lambda, Ta))
    {
        ioSynchronized_ = false;
        return;
    }
#endif

    // Prepare threads and scratch
    const bool parallel = useParallelCellLoops();
    const label nThreads = parallel ? maxCellLoopThreads() : 1;
    prepareSolveScratch(nThreads);

    // Execute the batched dispatch
    BatchedTensionBackend backend(*this);
    BatchedTensionExecutor<BatchedTensionBackend, CellScratch> executor
    (
        backend,
        parallel,
        solveScratch_
    );

    executor.solveCells(t, dt, driveSignals_, lambda, Ta);

    ioSynchronized_ = false;
}

} // End namespace Foam
