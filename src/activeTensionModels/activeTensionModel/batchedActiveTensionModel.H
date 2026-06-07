/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

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
        bool useRushLarsen_;
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

        virtual bool rushLarsenParametersForCell
        (
            const label cellI,
            const label stateI,
            const scalarUList& stateValues,
            const scalarUList& rateValues,
            const scalarUList& algebraicValues,
            scalar& steadyState,
            scalar& tau
        ) const
        {
            return false;
        }

        virtual scalar activeTensionFromHotPath
        (
            const label cellI,
            const scalarUList& stateValues,
            const scalarUList& algebraicValues
        ) const = 0;

        virtual scalar timeScaleFactor() const
        {
            return 1.0;
        }

        // Backend Bridge API for Executor
        friend class BatchedTensionBackend;
        struct BatchedTensionBackend
        {
            const batchedActiveTensionModel& model_;

            BatchedTensionBackend(const batchedActiveTensionModel& model) : model_(model) {}

            label nCells() const { return model_.nCells_; }
            label nSubsteps() const { return model_.nSubsteps_; }
            bool useRushLarsen() const { return model_.useRushLarsen_; }
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

            void applyExplicitEulerStep(const scalar dtSubstep, CellScratch& scratch) const
            {
                forAll(scratch.stateValues, i)
                {
                    scratch.stateValues[i] += dtSubstep * scratch.rateValues[i];
                }
            }

            void buildPredictorState(const label cellI, const scalar dtSubstep, CellScratch& scratch) const
            {
                forAll(scratch.stateValues, stateI)
                {
                    scalar steadyState = 0.0;
                    scalar tau = 0.0;
                    bool useRL = model_.rushLarsenParametersForCell
                    (
                        cellI, stateI, scratch.stateValues, scratch.rateValues,
                        scratch.algebraicValues, steadyState, tau
                    );

                    scratch.predictorRLApplies[stateI] = useRL ? 1.0 : 0.0;
                    if (useRL)
                    {
                        scratch.predictorStateValues[stateI] = steadyState + (scratch.stateValues[stateI] - steadyState) * std::exp(-dtSubstep/tau);
                    }
                    else
                    {
                        scratch.predictorStateValues[stateI] = scratch.stateValues[stateI] + dtSubstep * scratch.rateValues[stateI];
                    }
                }
            }

            void applyCorrectorState(CellScratch& scratch) const
            {
                forAll(scratch.stateValues, i)
                {
                    scratch.stateValues[i] = scratch.predictorStateValues[i];
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

    virtual ~batchedActiveTensionModel() = default;

    virtual void calculateTension
    (
        const scalar t,
        const scalar dt,
        const scalarField& lambda,
        scalarField& Ta
    ) override;

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
    virtual const PtrList<scalarField>* ioAlgebraicPtr() const override { syncAllToIO(); return &ioAlgebraics_; }
};

} // End namespace Foam

#endif
