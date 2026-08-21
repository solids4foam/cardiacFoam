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
    persistAlgebraics_(dict.lookupOrDefault<Switch>("storeBatchedAlgebraics", false)),
    useRushLarsen_(dict.lookupOrDefault<word>("batchedIntegrator", "euler") == "rushLarsen"),
    parallelCellUpdates_(dict.lookupOrDefault<Switch>("batchedParallelCells", false)),
    parallelMinCells_(dict.lookupOrDefault<label>("batchedParallelMinCells", 256)),
    core_(nIntegrationPoints, nStates, nAlgebraics),
    ioStates_(nIntegrationPoints),
    ioRates_(nIntegrationPoints),
    ioAlgebraics_(nIntegrationPoints),
    ioSynchronized_(false),
    solveScratch_(),
    solveScratchThreads_(-1)
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
            new CellScratch(nStates_, nAlgebraics_, useRushLarsen_)
        );
    }

    solveScratchThreads_ = nThreads;
}


void batchedActiveTensionModel::syncAllToIO() const
{
    if (ioSynchronized_) return;

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


void batchedActiveTensionModel::writeRestartState(const fvMesh& mesh) const
{
    syncAllToIO();
    restartStateIO::writeStates(mesh, type(), nStates_, ioStates_);
}


void batchedActiveTensionModel::refreshRestartState(const fvMesh& mesh)
{
    CellScratch scratch(nStates_, nAlgebraics_, useRushLarsen_);
    BatchedTensionBackend backend(*this);
    const ElectromechanicalSignalProvider& p = provider();

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        backend.gatherCellState(cellI, scratch.stateValues);
        scratch.resetPrimary();
        backend.evaluateScratchAtTime
        (
            cellI,
            mesh.time().value(),
            p.signal(cellI, driveSignal()),
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

    CellScratch scratch(nStates_, nAlgebraics_, useRushLarsen_);
    BatchedTensionBackend backend(*this);
    const ElectromechanicalSignalProvider& p = provider();

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        backend.gatherCellState(cellI, scratch.stateValues);
        scratch.resetPrimary();
        backend.evaluateScratchAtTime
        (
            cellI,
            0.0,
            p.signal(cellI, driveSignal()),
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

    const ElectromechanicalSignalProvider& p = provider();
    const CouplingSignal sig = driveSignal();

    if (driveSignals_.size() != nCells_)
    {
        driveSignals_.setSize(nCells_);
    }

    // Gather per-cell drive signals safely before parallel/GPU dispatch
    // This avoids calling virtual functions inside OpenMP/CUDA kernels
    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        driveSignals_[cellI] = p.signal(cellI, sig);
    }

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
