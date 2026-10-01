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

#ifndef batchedActiveTensionModel_H
#define batchedActiveTensionModel_H

#include "activeTensionModel.H"
#include "batchedActiveTensionCore.H"
#include "batchedActiveTensionExecution.H"
#ifdef HAS_CUDA
#include "batchedActiveTensionCoreCUDA.H"
#endif

namespace Foam
{

class batchedActiveTensionModel
:
    public activeTensionModel
{
protected:

    // Protected Data

        label nCells_;
        label nStates_;
        label nAlgebraics_;
        label nSubsteps_;
        bool persistAlgebraics_;
        Switch parallelCellUpdates_;
        label parallelMinCells_;

        mutable batchedActiveTensionCore core_;

        // IO storage pointers requested by the base class
        mutable PtrList<scalarField> ioStates_;
        mutable PtrList<scalarField> ioRates_;
        mutable PtrList<scalarField> ioAlgebraics_;
        mutable bool ioSynchronized_;

        using SolveFieldBuffer = scalarField;
        using CellScratch = BatchedTensionScratchT<SolveFieldBuffer>;

        mutable PtrList<CellScratch> solveScratch_;
        mutable label solveScratchThreads_;

        // Buffer to gather per-cell drive signals (Vm or Cai) before batch dispatch
        mutable scalarField driveSignals_;

#ifdef HAS_CUDA
        // Device execution is optional at runtime: a CUDA-built library still
        // falls back to the host path when no usable device is visible.
        mutable batchedActiveTensionCoreCUDA cuda_;
        Switch useCUDA_;
        Switch cudaProfile_;
        bool cudaEnabled_;
        mutable bool cudaAnnounced_;
        mutable bool cudaHostStateStale_;
        mutable label cudaProfileCalls_;
        mutable double cudaProfileH2DSeconds_;
        mutable double cudaProfileKernelSeconds_;
        mutable double cudaProfileD2HSeconds_;
        mutable double cudaProfileTotalSeconds_;
#endif

        // Private Member Functions

        static label maxCellLoopThreads()
        {
#ifdef _OPENMP
            return omp_get_max_threads();
#else
            return 1;
#endif
        }

        bool useParallelCellLoops() const
        {
            return parallelCellUpdates_ && nCells_ >= parallelMinCells_;
        }

        void prepareSolveScratch(const label nThreads) const;

        void syncAllToIO() const;
        void syncStatesFromIO() const;

        batchedActiveTensionModel(const batchedActiveTensionModel&) = delete;
        void operator=(const batchedActiveTensionModel&) = delete;

protected:

        batchedActiveTensionModel
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const label nStates,
            const label nAlgebraics
        );

        // Core Access
        inline scalar& state(const label cellI, const label stateI) const
        {
            return core_.stateRef(cellI, stateI);
        }

        inline scalar& rate(const label cellI, const label stateI) const
        {
            return core_.rateRef(cellI, stateI);
        }

        inline scalar& algebraic(const label cellI, const label algI) const
        {
            return core_.algebraicRef(cellI, algI);
        }

        // Implementation Hooks
        virtual void evaluateHotPathStateForCell
        (
            const label cellI,
            const scalar modelTime,
            const scalar driveSignal,
            const scalar lambda,
            const scalarUList& stateValues,
            scalarUList& rateValues,
            scalarUList& algebraicValues
        ) const = 0;

        virtual scalar activeTensionFromHotPath
        (
            const label cellI,
            const scalarUList& stateValues,
            const scalarUList& algebraicValues
        ) const = 0;

#ifdef HAS_CUDA
        // Fill model-specific input algebraics in the host SoA buffer before
        // it is mirrored to the device.  Land models require Cai, lambda and
        // lambda-rate; NashPanfilov takes its drive signal separately.
        virtual void prepareCUDAInputs
        (
            const scalarField& driveSignals,
            const scalarField& lambda
        ) const = 0;

        // Upload only time-varying input rows after the persistent state has
        // been established on the device.
        virtual void uploadCUDAInputs() const = 0;

        // Evaluate one model-specific right-hand side over all device cells.
        // The shared base applies the Euler state update between invocations.
        virtual void launchCUDAKernel() const = 0;

        // `Ta` is the only active-tension quantity required by the coupled
        // solver every step.  Download its model-specific state/algebraic row
        // without forcing all ODE state and diagnostic arrays back to host.
        virtual void downloadCUDATension(scalarField& Ta) const = 0;

        bool solveOnCUDA
        (
            const scalar t,
            const scalar dt,
            const scalarField& lambda,
            scalarField& Ta
        ) const;
#endif

        // timeScaleFactor() is declared pure virtual on activeTensionModel;
        // every concrete model states its own scale.

        // Backend Bridge API for Executor
        friend class BatchedTensionBackend;
        struct BatchedTensionBackend
        {
            const batchedActiveTensionModel& model_;

            BatchedTensionBackend(const batchedActiveTensionModel& model) : model_(model) {}

            label nCells() const { return model_.nCells_; }
            label nSubsteps() const { return model_.nSubsteps_; }
            scalar timeScaleFactor() const { return model_.timeScaleFactor(); }

            void gatherCellState(const label cellI, scalarField& stateValues) const
            {
                model_.core_.gatherCellState(cellI, stateValues);
            }

            void scatterCellState(const label cellI, const scalarField& stateValues) const
            {
                model_.core_.scatterCellState(cellI, stateValues);
            }

            void evaluateScratchAtTime
            (
                const label cellI,
                const scalar modelTime,
                const scalar driveSignal,
                const scalar lambda,
                CellScratch& scratch
            ) const
            {
                model_.evaluateHotPathStateForCell
                (
                    cellI, modelTime, driveSignal, lambda,
                    scratch.stateValues, scratch.rateValues, scratch.algebraicValues
                );
            }

            void advanceSubstep(const label cellI, const scalar dtSubstep, CellScratch& scratch) const
            {
                // Explicit Euler. Active-tension states are not stiff: the
                // fastest Land time constant is ~0.7 ms against a 0.01 ms
                // timestep, so the Rush-Larsen machinery the ionic models
                // need buys nothing here. batchedIntegrator is deliberately
                // not read: it is an ionic-model key, and singleCellSolver
                // hands both families the same dictionary.
                forAll(scratch.stateValues, stateI)
                {
                    scratch.stateValues[stateI] +=
                        dtSubstep*scratch.rateValues[stateI];
                }
            }

            void syncEvaluatedOutputs(const label cellI, const CellScratch& scratch) const
            {
                model_.core_.scatterCellRatesAndAlgebraics
                (
                    cellI, scratch.rateValues, scratch.algebraicValues, model_.persistAlgebraics_
                );
            }

            scalar activeTensionFromScratch(const label cellI, const CellScratch& scratch) const
            {
                return model_.activeTensionFromHotPath(cellI, scratch.stateValues, scratch.algebraicValues);
            }
        };

public:

    virtual ~batchedActiveTensionModel();

    virtual void calculateTension
    (
        const scalar t,
        const scalar dt,
        const scalarField& lambda,
        scalarField& Ta
    ) override;

    virtual bool supportsRestartState() const override;

    virtual bool readRestartState(const fvMesh& mesh) override;

    virtual void writeRestartState(const fvMesh& mesh) const override;

    virtual void refreshRestartState(const fvMesh& mesh) override;

    virtual bool setStates(const UList<scalarField>& states) override;

    virtual bool restartTension(scalarField& Ta) const override;

    // Optional legacy override, but calculateTension is the real entry point
    virtual void solveAtPoint
    (
        const label integrationPtI,
        const scalar driveSignal,
        const scalar lambda,
        scalar& Ta
    ) override
    {
        // Not used by calculateTension in batched models
    }

    // Scalar derivatives override (unused by batched architecture)
    virtual void derivatives
    (
        const scalar t,
        const scalarField& y,
        scalarField& dydt
    ) const override
    {
        FatalErrorInFunction
            << "Batched models do not use the scalar derivatives() function. "
            << "They execute via evaluateHotPathStateForCell."
            << abort(FatalError);
    }

    // IO Support Hooks
    virtual label ioNumStates() const override { return nStates_; }
    virtual label ioNumAlgebraic() const override { return nAlgebraics_; }
    virtual const PtrList<scalarField>* ioStatesPtr() const override { syncAllToIO(); return &ioStates_; }
    virtual const PtrList<scalarField>* ioRatesPtr() const override { syncAllToIO(); return &ioRates_; }
    virtual const PtrList<scalarField>* ioAlgebraicPtr() const override
    {
        if (!persistAlgebraics_)
        {
            return nullptr;
        }
        syncAllToIO();
        return &ioAlgebraics_;
    }
};

} // End namespace Foam

#endif
