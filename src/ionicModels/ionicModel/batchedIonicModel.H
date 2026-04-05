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
    Foam::batchedIonicModel

Description
    Reusable base for GPU-friendly ionic models that execute as batched,
    explicit cell updates while preserving the existing ionicModel API.
\*---------------------------------------------------------------------------*/

#ifndef batchedIonicModel_H
#define batchedIonicModel_H

#include "ionicModelGPU.H"
#include "../ionicModel/batchedIonicCore.H"
#include "../ionicModel/batchedKernelExecution.H"
#include "PstreamReduceOps.H"
#include <chrono>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Foam
{

class batchedIonicModel
:
    public ionicModelGPU
{
    template<class Model>
    friend class BatchedKernelModelHooksView;
    template<class StorageBackend, class ModelPolicy, class FieldBuffer, class StepBuffer>
    friend class BatchedKernelHostBackend;

    // Private Data

        label nCells_;
        label nStates_;
        label nAlgebraics_;
        label nHotPathSupport_;
        label nSubsteps_;
        word integratorName_;
        Switch profileBatchedModel_;
        Switch parallelCellUpdates_;
        label parallelMinCells_;
        bool persistAlgebraics_;

        mutable batchedIonicCore core_;

        mutable PtrList<scalarField> ioStates_;
        mutable PtrList<scalarField> ioRates_;
        mutable PtrList<scalarField> ioAlgebraics_;
        mutable bool ioSynchronized_;
        mutable bool transientAlgebraicsSynchronized_;
        mutable scalar ioEvaluationModelTime_;
        mutable scalar totalSolveODETime_;
        mutable scalar totalVmInputSyncTime_;
        mutable scalar totalGatherScatterTime_;
        mutable scalar totalEvaluateStateTime_;
        mutable scalar totalStateUpdateTime_;
        mutable scalar totalStabilizeTime_;
        mutable scalar totalFinalSyncTime_;
        mutable label totalEvaluateStateCalls_;
        mutable label totalSolveODECalls_;


    // Private Member Functions

        using SolveFieldBuffer = scalarField;
        using SolveStepBuffer = scalarList;
        using SolveProfileTotals = BatchedKernelProfileTotals;
        using SolveKernelConfig = BatchedKernelConfig;
        using SolveStorageBackend =
            BatchedKernelSoAStorageView<SolveFieldBuffer, SolveStepBuffer>;
        using CellScratch = BatchedKernelScratchT<SolveFieldBuffer>;
        using SolveModelHooks = BatchedKernelModelHooksView<batchedIonicModel>;

        using SolveModelPolicy =
            BatchedKernelModelPolicy<SolveModelHooks, SolveFieldBuffer>;

        mutable PtrList<CellScratch> solveScratch_;
        mutable label solveScratchThreads_;
        mutable label solveScratchSupportSize_;
        mutable bool solveScratchPredictorState_;
        mutable bool solveScratchPredictorEvaluation_;
        mutable bool solveScratchSupportOnlyHotPath_;

        static label maxCellLoopThreads()
        {
#ifdef _OPENMP
            return omp_get_max_threads();
#else
            return 1;
#endif
        }

        static label cellLoopThreadIndex()
        {
#ifdef _OPENMP
            return omp_get_thread_num();
#else
            return 0;
#endif
        }

        bool useParallelCellLoops() const
        {
            return parallelCellUpdates_ && nCells_ >= parallelMinCells_;
        }

        batchedIonicModel(const batchedIonicModel&) = delete;
        void operator=(const batchedIonicModel&) = delete;

        bool persistentAlgebraicsRequested() const
        {
            return dict().lookupOrDefault<Switch>
            (
                "storeBatchedAlgebraics",
                false
            );
        }


protected:

        batchedIonicModel
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver,
            const label nStates,
            const label nAlgebraics
        )
        :
            ionicModelGPU
            (
                dict,
                nIntegrationPoints,
                initialDeltaT,
                solveVmWithinODESolver
            ),
            nCells_(nIntegrationPoints),
            nStates_(nStates),
            nAlgebraics_(nAlgebraics),
            nHotPathSupport_(0),
            nSubsteps_(dict.lookupOrDefault<label>("batchedSubsteps", 1)),
            integratorName_
            (
                dict.lookupOrDefault<word>("batchedIntegrator", "rushLarsenHeun")
            ),
            profileBatchedModel_
            (
                dict.lookupOrDefault<Switch>("profileBatchedModel", false)
            ),
            parallelCellUpdates_
            (
                dict.lookupOrDefault<Switch>("batchedParallelCells", false)
            ),
            parallelMinCells_
            (
                dict.lookupOrDefault<label>("batchedParallelMinCells", 256)
            ),
            persistAlgebraics_(true),
            core_(nIntegrationPoints, nStates, nAlgebraics),
            ioStates_(nIntegrationPoints),
            ioRates_(nIntegrationPoints),
            ioAlgebraics_(nIntegrationPoints),
            ioSynchronized_(false),
            transientAlgebraicsSynchronized_(false),
            ioEvaluationModelTime_(0.0),
            totalSolveODETime_(0.0),
            totalVmInputSyncTime_(0.0),
            totalGatherScatterTime_(0.0),
            totalEvaluateStateTime_(0.0),
            totalStateUpdateTime_(0.0),
            totalStabilizeTime_(0.0),
            totalFinalSyncTime_(0.0),
            totalEvaluateStateCalls_(0),
            totalSolveODECalls_(0),
            solveScratch_(),
            solveScratchThreads_(-1),
            solveScratchSupportSize_(-1),
            solveScratchPredictorState_(false),
            solveScratchPredictorEvaluation_(false),
            solveScratchSupportOnlyHotPath_(false)
        {
            if (nSubsteps_ < 1)
            {
                nSubsteps_ = 1;
            }

            if (parallelMinCells_ < 1)
            {
                parallelMinCells_ = 1;
            }

#ifndef _OPENMP
            if (parallelCellUpdates_)
            {
                WarningInFunction
                    << "batchedParallelCells was requested for ionic model "
                    << type()
                    << " but this library was built without OpenMP support. "
                    << "Falling back to serial cell loops." << nl;
                parallelCellUpdates_ = false;
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

        void configurePersistentAlgebraics()
        {
            const bool requested =
                dict().lookupOrDefault<Switch>
                (
                    "storeBatchedAlgebraics",
                    persistentAlgebraicsRequested()
                );

            if (requested == persistAlgebraics_)
            {
                return;
            }

            persistAlgebraics_ = requested;
            ioSynchronized_ = false;

            if (persistAlgebraics_)
            {
                core_.setAlgebraicsStorageEnabled(true);
                ioAlgebraics_.setSize(nCells_);

                for (label cellI = 0; cellI < nCells_; ++cellI)
                {
                    ioAlgebraics_.set(cellI, new scalarField(nAlgebraics_, 0.0));
                }
            }
            else
            {
                core_.setAlgebraicsStorageEnabled(false);
                ioAlgebraics_.clear();
            }

            transientAlgebraicsSynchronized_ = false;
        }

        void prepareSolveScratch
        (
            const label nThreads,
            const label supportSize,
            const bool usePredictorState,
            const bool usePredictorEvaluation,
            const bool useSupportOnlyHotPath
        ) const
        {
            if
            (
                solveScratchThreads_ == nThreads
             && solveScratchSupportSize_ == supportSize
             && solveScratchPredictorState_ == usePredictorState
             && solveScratchPredictorEvaluation_ == usePredictorEvaluation
             && solveScratchSupportOnlyHotPath_ == useSupportOnlyHotPath
             && solveScratch_.size() == nThreads
            )
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
                    new CellScratch
                    (
                        nStates_,
                        nAlgebraics_,
                        supportSize,
                        usePredictorState,
                        usePredictorEvaluation,
                        useSupportOnlyHotPath
                    )
                );
            }

            solveScratchThreads_ = nThreads;
            solveScratchSupportSize_ = supportSize;
            solveScratchPredictorState_ = usePredictorState;
            solveScratchPredictorEvaluation_ = usePredictorEvaluation;
            solveScratchSupportOnlyHotPath_ = useSupportOnlyHotPath;
        }

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

        inline scalar& support(const label cellI, const label supportI) const
        {
            return core_.supportRef(cellI, supportI);
        }

        inline label nCells() const
        {
            return nCells_;
        }

        inline label nSubsteps() const
        {
            return nSubsteps_;
        }

        inline label nHotPathSupport() const
        {
            return nHotPathSupport_;
        }

        void setHotPathSupportSize(const label nSupport)
        {
            nHotPathSupport_ = max(nSupport, label(0));
            core_.setHotPathSupportSize(nHotPathSupport_);
            solveScratchSupportSize_ = -1;
        }

        void gatherCellState
        (
            const label cellI,
            scalarField& stateValues
        ) const
        {
            core_.gatherCellState(cellI, stateValues);
        }

        void scatterCellState
        (
            const label cellI,
            const scalarField& stateValues
        ) const
        {
            core_.scatterCellState(cellI, stateValues);
        }

        void scatterCellRatesAndAlgebraics
        (
            const label cellI,
            const scalarUList& rateValues,
            const scalarUList& algebraicValues
        ) const
        {
            core_.scatterCellRatesAndAlgebraics
            (
                cellI,
                rateValues,
                algebraicValues,
                persistAlgebraics_
            );
        }

        void scatterCellRates
        (
            const label cellI,
            const scalarUList& rateValues
        ) const
        {
            core_.scatterCellRates(cellI, rateValues);
        }

        void scatterCellSupport
        (
            const label cellI,
            const scalarUList& supportValues
        ) const
        {
            core_.scatterCellSupport(cellI, supportValues);
        }

        void syncCellToIO(const label cellI) const
        {
            scalarField& stateIO = ioStates_[cellI];
            scalarField& rateIO = ioRates_[cellI];

            for (label stateI = 0; stateI < nStates_; ++stateI)
            {
                stateIO[stateI] = state(cellI, stateI);
                rateIO[stateI] = rate(cellI, stateI);
            }

            if (!persistAlgebraics_)
            {
                return;
            }

            scalarField& algIO = ioAlgebraics_[cellI];
            for (label algI = 0; algI < nAlgebraics_; ++algI)
            {
                algIO[algI] = algebraic(cellI, algI);
            }
        }

        void syncAllToIO() const
        {
            for (label cellI = 0; cellI < nCells_; ++cellI)
            {
                syncCellToIO(cellI);
            }

            ioSynchronized_ = true;
            transientAlgebraicsSynchronized_ = persistAlgebraics_;
        }

        void markIODirty() const
        {
            ioSynchronized_ = false;
            transientAlgebraicsSynchronized_ = false;
        }

        void setIOEvaluationModelTime(const scalar modelTime) const
        {
            ioEvaluationModelTime_ = modelTime;
        }

        void ensureTransientAlgebraicsSynchronized() const
        {
            if (transientAlgebraicsSynchronized_)
            {
                return;
            }

            ensureIOSynchronized();

            if (ioAlgebraics_.size() != nCells_)
            {
                ioAlgebraics_.setSize(nCells_);

                for (label cellI = 0; cellI < nCells_; ++cellI)
                {
                    ioAlgebraics_.set(cellI, new scalarField(nAlgebraics_, 0.0));
                }
            }

            scalarField stateValues(nStates_, 0.0);
            scalarField rateValues(nStates_, 0.0);
            scalarField algebraicValues(nAlgebraics_, 0.0);

            for (label cellI = 0; cellI < nCells_; ++cellI)
            {
                gatherCellState(cellI, stateValues);
                rateValues = 0.0;
                algebraicValues = 0.0;

                evaluateState
                (
                    ioEvaluationModelTime_,
                    stateValues,
                    rateValues,
                    algebraicValues
                );

                ioAlgebraics_[cellI] = algebraicValues;
            }

            transientAlgebraicsSynchronized_ = true;
        }

        void ensureIOSynchronized() const
        {
            if (!ioSynchronized_)
            {
                syncAllToIO();
            }
        }

        void syncStatesFromIO() const
        {
            for (label cellI = 0; cellI < nCells_; ++cellI)
            {
                const scalarField& stateIO = ioStates_[cellI];
                const scalarField& rateIO = ioRates_[cellI];

                for (label stateI = 0; stateI < nStates_; ++stateI)
                {
                    state(cellI, stateI) = stateIO[stateI];
                    rate(cellI, stateI) = rateIO[stateI];
                }

                if (!persistAlgebraics_)
                {
                    continue;
                }

                const scalarField& algIO = ioAlgebraics_[cellI];
                for (label algI = 0; algI < nAlgebraics_; ++algI)
                {
                    algebraic(cellI, algI) = algIO[algI];
                }
            }
        }

        void clearTransientIOData() const
        {
            core_.clearTransientSolveData(persistAlgebraics_);
            transientAlgebraicsSynchronized_ = false;

            if (!persistAlgebraics_)
            {
                ioAlgebraics_.clear();
            }

            for (label cellI = 0; cellI < nCells_; ++cellI)
            {
                ioRates_[cellI] = 0.0;

                if (persistAlgebraics_)
                {
                    ioAlgebraics_[cellI] = 0.0;
                }
            }
        }

        void evaluateScratchAtTime
        (
            const scalar modelTime,
            CellScratch& scratch,
            const SolveKernelConfig& kernel
        ) const
        {
            if (kernel.useSupportOnlyHotPath)
            {
                evaluateHotPathState
                (
                    modelTime,
                    scratch.stateValues,
                    scratch.rateValues,
                    scratch.supportValues
                );
            }
            else
            {
                evaluateState
                (
                    modelTime,
                    scratch.stateValues,
                    scratch.rateValues,
                    scratch.algebraicValues
                );
            }

            if (!kernel.solveVm)
            {
                scratch.rateValues[kernel.vmStateI] = 0.0;
            }
        }

        void evaluatePredictorScratchAtTime
        (
            const scalar modelTime,
            CellScratch& scratch,
            const SolveKernelConfig& kernel
        ) const
        {
            if (kernel.useSupportOnlyHotPath)
            {
                evaluateHotPathState
                (
                    modelTime,
                    scratch.predictorStateValues,
                    scratch.predictorRateValues,
                    scratch.predictorSupportValues
                );
            }
            else
            {
                evaluateState
                (
                    modelTime,
                    scratch.predictorStateValues,
                    scratch.predictorRateValues,
                    scratch.predictorAlgebraicValues
                );
            }

            if (!kernel.solveVm)
            {
                scratch.predictorRateValues[kernel.vmStateI] = 0.0;
            }
        }

        bool rushLarsenEnabledForState
        (
            const label stateI,
            const CellScratch& scratch,
            const SolveKernelConfig& kernel,
            scalar& steadyState,
            scalar& tau
        ) const
        {
            if (!kernel.useRushLarsen)
            {
                return false;
            }

            const bool hasParameters =
                kernel.useSupportOnlyHotPath
              ? rushLarsenParametersFromHotPathSupport
                (
                    stateI,
                    scratch.stateValues,
                    scratch.rateValues,
                    scratch.supportValues,
                    steadyState,
                    tau
                )
              : rushLarsenParameters
                (
                    stateI,
                    scratch.stateValues,
                    scratch.rateValues,
                    scratch.algebraicValues,
                    steadyState,
                    tau
                );

            return hasParameters
                && tau > VSMALL
                && std::isfinite(steadyState)
                && std::isfinite(tau);
        }

        void buildPredictorState
        (
            CellScratch& scratch,
            const SolveKernelConfig& kernel
        ) const
        {
            for (label stateI = 0; stateI < nStates_; ++stateI)
            {
                if (!kernel.solveVm && stateI == kernel.vmStateI)
                {
                    scratch.predictorStateValues[stateI] =
                        scratch.stateValues[stateI];
                    continue;
                }

                scalar steadyState = 0.0;
                scalar tau = 0.0;
                const bool useRLForState =
                    rushLarsenEnabledForState
                    (
                        stateI,
                        scratch,
                        kernel,
                        steadyState,
                        tau
                    );

                if (useRLForState)
                {
                    scratch.predictorStateValues[stateI] =
                        steadyState
                      + (scratch.stateValues[stateI] - steadyState)
                       * std::exp(-kernel.dtSubstep/tau);
                }
                else
                {
                    scratch.predictorStateValues[stateI] =
                        scratch.stateValues[stateI]
                      + kernel.dtSubstep*scratch.rateValues[stateI];
                }
            }
        }

        void applyCorrectorState
        (
            CellScratch& scratch,
            const SolveKernelConfig& kernel
        ) const
        {
            for (label stateI = 0; stateI < nStates_; ++stateI)
            {
                if (!kernel.solveVm && stateI == kernel.vmStateI)
                {
                    continue;
                }

                scalar steadyState = 0.0;
                scalar tau = 0.0;
                const bool useRLForState =
                    rushLarsenEnabledForState
                    (
                        stateI,
                        scratch,
                        kernel,
                        steadyState,
                        tau
                    );

                if (useRLForState)
                {
                    scratch.stateValues[stateI] =
                        scratch.predictorStateValues[stateI];
                }
                else if (kernel.useHeun)
                {
                    scratch.stateValues[stateI] +=
                        0.5*kernel.dtSubstep
                       *
                        (
                            scratch.rateValues[stateI]
                          + scratch.predictorRateValues[stateI]
                        );
                }
                else
                {
                    scratch.stateValues[stateI] =
                        scratch.predictorStateValues[stateI];
                }
            }
        }

        void syncEvaluatedOutputs
        (
            const label cellI,
            const CellScratch& scratch,
            const SolveKernelConfig& kernel
        ) const
        {
            if (kernel.useSupportOnlyHotPath)
            {
                scatterCellRates(cellI, scratch.rateValues);
                scatterCellSupport(cellI, scratch.supportValues);
            }
            else
            {
                scatterCellRatesAndAlgebraics
                (
                    cellI,
                    scratch.rateValues,
                    scratch.algebraicValues
                );
            }
        }

        scalar ionicCurrentFromScratch
        (
            const CellScratch& scratch,
            const SolveKernelConfig& kernel
        ) const
        {
            if (kernel.useSupportOnlyHotPath)
            {
                return ionicCurrentFromHotPathSupport(scratch.supportValues);
            }

            return ionicCurrentFromEvaluation(scratch.algebraicValues);
        }

        bool persistentAlgebraicsEnabled() const
        {
            return persistAlgebraics_;
        }

        virtual label activeHotPathSupportSize
        (
            const bool useExplicitEuler
        ) const
        {
            return nHotPathSupport_;
        }

        virtual label voltageStateIndex() const override = 0;
        virtual scalar timeScaleFactor() const = 0;
        virtual scalar vmToState(const scalar vm) const = 0;
        virtual scalar ionicCurrentFromEvaluation
        (
            const scalarUList& algebraicValues
        ) const = 0;

        virtual scalar ionicCurrentFromHotPathSupport
        (
            const scalarUList& supportValues
        ) const
        {
            return ionicCurrentFromEvaluation(supportValues);
        }

        virtual bool usesCompactSupport() const
        {
            return false;
        }

        virtual void evaluateState
        (
            const scalar modelTime,
            const scalarUList& stateValues,
            scalarUList& rateValues,
            scalarUList& algebraicValues
        ) const = 0;

        virtual void evaluateHotPathState
        (
            const scalar modelTime,
            const scalarUList& stateValues,
            scalarUList& rateValues,
            scalarUList& supportValues
        ) const
        {
            evaluateState(modelTime, stateValues, rateValues, supportValues);
        }

        virtual bool rushLarsenParameters
        (
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

        virtual bool rushLarsenParametersFromHotPathSupport
        (
            const label stateI,
            const scalarUList& stateValues,
            const scalarUList& rateValues,
            const scalarUList& supportValues,
            scalar& steadyState,
            scalar& tau
        ) const
        {
            return rushLarsenParameters
            (
                stateI,
                stateValues,
                rateValues,
                supportValues,
                steadyState,
                tau
            );
        }

        virtual void stabilizeStateValues(scalarUList& stateValues) const
        {}

        virtual void stabilizeCellState(const label cellI) const
        {}

        virtual const PtrList<scalarField>* ioStatesPtr() const override
        {
            ensureIOSynchronized();
            return &ioStates_;
        }

        virtual const PtrList<scalarField>* ioAlgebraicPtr() const override
        {
            if (!transientAlgebraicsSynchronized_)
            {
                return nullptr;
            }

            ensureIOSynchronized();
            return &ioAlgebraics_;
        }

        virtual void prepareIOAccess
        (
            const wordList&,
            const bool needsAlgebraics
        ) const override
        {
            if (needsAlgebraics)
            {
                ensureTransientAlgebraicsSynchronized();
            }
        }

        virtual const PtrList<scalarField>* ioRatesPtr() const override
        {
            ensureIOSynchronized();
            return &ioRates_;
        }

        virtual label ioNumStates() const override
        {
            return nStates_;
        }

        virtual label ioNumAlgebraic() const override
        {
            return nAlgebraics_;
        }

        virtual label gpuCellCount() const override
        {
            return nCells_;
        }

        virtual scalar gpuStateValue
        (
            const label cellI,
            const label stateI
        ) const override
        {
            return state(cellI, stateI);
        }

        virtual scalar gpuRateValue
        (
            const label cellI,
            const label stateI
        ) const override
        {
            return rate(cellI, stateI);
        }

        virtual void gpuGatherCellState
        (
            const label cellI,
            scalarField& stateValues
        ) const override
        {
            gatherCellState(cellI, stateValues);
        }

        template<class Model>
        void solveODEImpl
        (
            const Model& model,
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im
        )
        {
            using ModelHooks = BatchedKernelModelHooksView<Model>;
            using ModelPolicy =
                BatchedKernelModelPolicy<ModelHooks, SolveFieldBuffer>;
            using clock = std::chrono::steady_clock;

            const bool profileSolve = profileBatchedModel_;
            const bool emitDebugSample = debugPrintedNames().size();
            const auto solveStart =
                profileSolve ? clock::now() : clock::time_point();
            const scalar tStart = stepStartTime*timeScaleFactor();
            const scalar dtModel = deltaT*timeScaleFactor();
            const scalar tEnd = tStart + dtModel;
            const scalar dtSubstep = dtModel/scalar(nSubsteps_);
            const label vmStateI = voltageStateIndex();
            SolveKernelConfig kernel;
            kernel.dtSubstep = dtSubstep;
            kernel.vmStateI = vmStateI;
            kernel.solveVm = solveVmWithinODESolver();
            kernel.useHeun = integratorName_ == "heun"
                || integratorName_ == "rushLarsenHeun"
                || integratorName_ == "mixed";
            kernel.useRushLarsen = integratorName_ == "rushLarsenHeun"
                || integratorName_ == "mixed";
            kernel.useExplicitEuler =
                !kernel.useHeun && !kernel.useRushLarsen;
            const label sampleCell =
                emitDebugSample ? sampleIntegrationPoint(nCells_) : -1;
            kernel.useSupportOnlyHotPath =
                usesCompactSupport()
             && nHotPathSupport_ > 0
             && !utilitiesMode();
            const label hotPathSupportSize =
                kernel.useSupportOnlyHotPath
              ? activeHotPathSupportSize(kernel.useExplicitEuler)
              : nHotPathSupport_;
            const bool parallelCellLoops = useParallelCellLoops();
            const label nCellLoopThreads =
                parallelCellLoops ? maxCellLoopThreads() : 1;
            List<SolveProfileTotals> threadProfiles(nCellLoopThreads);

            prepareSolveScratch
            (
                nCellLoopThreads,
                hotPathSupportSize,
                !kernel.useExplicitEuler,
                kernel.useHeun,
                kernel.useSupportOnlyHotPath
            );

            SolveStorageBackend storageBackend
            (
                nCells_,
                nStates_,
                core_.statesField(),
                step()
            );
            ModelHooks modelHooks(model);
            ModelPolicy modelPolicy(modelHooks);
            BatchedKernelHostBackend
            <
                SolveStorageBackend,
                ModelPolicy,
                SolveFieldBuffer,
                SolveStepBuffer
            > backend(storageBackend, modelPolicy);
            BatchedKernelExecutor
            <
                BatchedKernelHostBackend
                <
                    SolveStorageBackend,
                    ModelPolicy,
                    SolveFieldBuffer,
                    SolveStepBuffer
                >,
                SolveKernelConfig,
                CellScratch,
                SolveProfileTotals
            > executor
            (
                backend,
                kernel,
                parallelCellLoops,
                profileSolve,
                solveScratch_,
                threadProfiles
            );

            executor.syncVmInputs(Vm);

            markIODirty();
            setIOEvaluationModelTime(tStart);

            executor.emitDebugSample(sampleCell, tStart, tEnd, dtModel);
            executor.advanceSubsteps(tStart, tEnd, Im);

            setIOEvaluationModelTime(tEnd);

            if (sampleCell >= 0 && sampleCell < nCells_)
            {
                debugPrintFields(sampleCell, tStart, tEnd, dtModel);
            }

            forAll(threadProfiles, threadI)
            {
                totalVmInputSyncTime_ += threadProfiles[threadI].vmInputSyncTime;
                totalGatherScatterTime_ += threadProfiles[threadI].gatherScatterTime;
                totalEvaluateStateTime_ += threadProfiles[threadI].evaluateStateTime;
                totalStateUpdateTime_ += threadProfiles[threadI].stateUpdateTime;
                totalStabilizeTime_ += threadProfiles[threadI].stabilizeTime;
                totalFinalSyncTime_ += threadProfiles[threadI].finalSyncTime;
                totalEvaluateStateCalls_ += threadProfiles[threadI].evaluateStateCalls;
            }

            if (profileSolve)
            {
                totalSolveODETime_ +=
                    std::chrono::duration<scalar>(clock::now() - solveStart).count();
                ++totalSolveODECalls_;
            }
        }

public:

        virtual void writeProfilingSummary() const override
        {
            if (!profileBatchedModel_)
            {
                return;
            }

            auto printScalarMetric = [&](const word& metricName, const scalar localValue)
            {
                scalar sumValue = localValue;
                scalar maxValue = localValue;
                reduce(sumValue, sumOp<scalar>());
                reduce(maxValue, maxOp<scalar>());

                if (Pstream::master())
                {
                    Info<< "PROFILE_BATCHED model=" << type()
                        << " metric=" << metricName
                        << " sum=" << sumValue
                        << " max=" << maxValue << nl;
                }
            };

            auto printLabelMetric = [&](const word& metricName, const label localValue)
            {
                label sumValue = localValue;
                label maxValue = localValue;
                reduce(sumValue, sumOp<label>());
                reduce(maxValue, maxOp<label>());

                if (Pstream::master())
                {
                    Info<< "PROFILE_BATCHED model=" << type()
                        << " metric=" << metricName
                        << " sum=" << sumValue
                        << " max=" << maxValue << nl;
                }
            };

            printScalarMetric("totalSolveODETime_s", totalSolveODETime_);
            printScalarMetric("vmInputSyncTime_s", totalVmInputSyncTime_);
            printScalarMetric("gatherScatterTime_s", totalGatherScatterTime_);
            printScalarMetric("evaluateStateTime_s", totalEvaluateStateTime_);
            printScalarMetric("stateUpdateTime_s", totalStateUpdateTime_);
            printScalarMetric("stabilizeTime_s", totalStabilizeTime_);
            printScalarMetric("finalSyncTime_s", totalFinalSyncTime_);
            printLabelMetric("evaluateStateCalls", totalEvaluateStateCalls_);
            printLabelMetric("solveODECalls", totalSolveODECalls_);
        }

        virtual label nEqns() const override
        {
            return nStates_;
        }

        virtual void importFields
        (
            const volScalarField& Vm,
            const wordList& fieldNames,
            const PtrList<volScalarField>& inFields
        ) override
        {
            ensureIOSynchronized();
            ionicModel::importFields(Vm, fieldNames, inFields);
            clearTransientIOData();
            syncStatesFromIO();
            ioSynchronized_ = true;
        }

        virtual void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im
        ) override
        {
            solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
        }
};

} // End namespace Foam

#endif

// ************************************************************************* //
