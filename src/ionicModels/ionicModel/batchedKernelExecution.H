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
    Foam::BatchedKernelHostBackend
    Foam::BatchedKernelExecutor

Description
    Reusable batched-kernel execution layer for SoA ionic stepping.
\*---------------------------------------------------------------------------*/

#ifndef batchedKernelExecution_H
#define batchedKernelExecution_H

#include <chrono>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Foam
{

struct BatchedKernelProfileTotals
{
    scalar vmInputSyncTime;
    scalar gatherScatterTime;
    scalar evaluateStateTime;
    scalar stateUpdateTime;
    scalar stabilizeTime;
    scalar finalSyncTime;
    label evaluateStateCalls;

    BatchedKernelProfileTotals()
    :
        vmInputSyncTime(0.0),
        gatherScatterTime(0.0),
        evaluateStateTime(0.0),
        stateUpdateTime(0.0),
        stabilizeTime(0.0),
        finalSyncTime(0.0),
        evaluateStateCalls(0)
    {}
};


struct BatchedKernelConfig
{
    scalar dtSubstep;
    label vmStateI;
    bool solveVm;
    bool useExplicitEuler;
    bool useHeun;
    bool useRushLarsen;
    bool useSupportOnlyHotPath;

    BatchedKernelConfig()
    :
        dtSubstep(0.0),
        vmStateI(-1),
        solveVm(false),
        useExplicitEuler(false),
        useHeun(false),
        useRushLarsen(false),
        useSupportOnlyHotPath(false)
    {}
};


template<class FieldBuffer, class StepBuffer>
class BatchedKernelSoAStorageView
{
    label nCells_;
    label nStates_;
    FieldBuffer& statesSoA_;
    StepBuffer& stepField_;

    inline label offset(const label cellI, const label stateI) const
    {
        return stateI*nCells_ + cellI;
    }

public:

    BatchedKernelSoAStorageView
    (
        const label nCells,
        const label nStates,
        FieldBuffer& statesSoA,
        StepBuffer& stepField
    )
    :
        nCells_(nCells),
        nStates_(nStates),
        statesSoA_(statesSoA),
        stepField_(stepField)
    {}

    label nCells() const
    {
        return nCells_;
    }

    scalar& stateRef(const label cellI, const label stateI) const
    {
        return statesSoA_[offset(cellI, stateI)];
    }

    void gatherCellState
    (
        const label cellI,
        FieldBuffer& stateValues
    ) const
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            stateValues[stateI] = statesSoA_[offset(cellI, stateI)];
        }
    }

    void scatterCellState
    (
        const label cellI,
        const FieldBuffer& stateValues
    ) const
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            statesSoA_[offset(cellI, stateI)] = stateValues[stateI];
        }
    }

    StepBuffer& stepField() const
    {
        return stepField_;
    }
};


template<class FieldBuffer>
struct BatchedKernelScratchT
{
    FieldBuffer stateValues;
    FieldBuffer rateValues;
    FieldBuffer algebraicValues;
    FieldBuffer supportValues;
    FieldBuffer predictorStateValues;
    FieldBuffer predictorRateValues;
    FieldBuffer predictorAlgebraicValues;
    FieldBuffer predictorSupportValues;

    BatchedKernelScratchT
    (
        const label nStates,
        const label nAlgebraics,
        const label nSupport,
        const bool usePredictorState,
        const bool usePredictorEvaluation,
        const bool useSupportOnlyHotPath
    )
    :
        stateValues(nStates, 0.0),
        rateValues(nStates, 0.0),
        algebraicValues(useSupportOnlyHotPath ? 0 : nAlgebraics, 0.0),
        supportValues(nSupport, 0.0),
        predictorStateValues(usePredictorState ? nStates : 0, 0.0),
        predictorRateValues(usePredictorEvaluation ? nStates : 0, 0.0),
        predictorAlgebraicValues
        (
            usePredictorEvaluation && !useSupportOnlyHotPath ? nAlgebraics : 0,
            0.0
        ),
        predictorSupportValues(usePredictorEvaluation ? nSupport : 0, 0.0)
    {}

    void resetPrimary(const bool useSupportOnlyHotPath)
    {
        rateValues = 0.0;

        if (useSupportOnlyHotPath)
        {
            supportValues = 0.0;
        }
        else
        {
            algebraicValues = 0.0;
        }
    }

    void resetPredictor
    (
        const bool usePredictorEvaluation,
        const bool useSupportOnlyHotPath
    )
    {
        if (!usePredictorEvaluation)
        {
            return;
        }

        predictorRateValues = 0.0;

        if (useSupportOnlyHotPath)
        {
            predictorSupportValues = 0.0;
        }
        else
        {
            predictorAlgebraicValues = 0.0;
        }
    }
};


template<class Model>
class BatchedKernelModelHooksView
{
    const Model& model_;

public:

    explicit BatchedKernelModelHooksView(const Model& model)
    :
        model_(model)
    {}

    label nSubsteps() const
    {
        return model_.nSubsteps_;
    }

    scalar vmInputToState(const scalar vm) const
    {
        return model_.vmToState(vm);
    }

    void evaluateState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& algebraicValues
    ) const
    {
        model_.evaluateState
        (
            modelTime,
            stateValues,
            rateValues,
            algebraicValues
        );
    }

    void evaluateHotPathState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& supportValues
    ) const
    {
        model_.evaluateHotPathState
        (
            modelTime,
            stateValues,
            rateValues,
            supportValues
        );
    }

    bool rushLarsenParameters
    (
        const label stateI,
        const scalarUList& stateValues,
        const scalarUList& rateValues,
        const scalarUList& algebraicValues,
        scalar& steadyState,
        scalar& tau
    ) const
    {
        return model_.rushLarsenParameters
        (
            stateI,
            stateValues,
            rateValues,
            algebraicValues,
            steadyState,
            tau
        );
    }

    bool rushLarsenParametersFromHotPathSupport
    (
        const label stateI,
        const scalarUList& stateValues,
        const scalarUList& rateValues,
        const scalarUList& supportValues,
        scalar& steadyState,
        scalar& tau
    ) const
    {
        return model_.rushLarsenParametersFromHotPathSupport
        (
            stateI,
            stateValues,
            rateValues,
            supportValues,
            steadyState,
            tau
        );
    }

    void scatterCellRates
    (
        const label cellI,
        const scalarUList& rateValues
    ) const
    {
        model_.scatterCellRates(cellI, rateValues);
    }

    void scatterCellSupport
    (
        const label cellI,
        const scalarUList& supportValues
    ) const
    {
        model_.scatterCellSupport(cellI, supportValues);
    }

    void scatterCellRatesAndAlgebraics
    (
        const label cellI,
        const scalarUList& rateValues,
        const scalarUList& algebraicValues
    ) const
    {
        model_.scatterCellRatesAndAlgebraics
        (
            cellI,
            rateValues,
            algebraicValues
        );
    }

    scalar ionicCurrentFromHotPathSupport
    (
        const scalarUList& supportValues
    ) const
    {
        return model_.ionicCurrentFromHotPathSupport(supportValues);
    }

    scalar ionicCurrentFromEvaluation
    (
        const scalarUList& algebraicValues
    ) const
    {
        return model_.ionicCurrentFromEvaluation(algebraicValues);
    }

    void stabilizeCellState(const label cellI) const
    {
        model_.stabilizeCellState(cellI);
    }

    void stabilizeStateValues(scalarUList& stateValues) const
    {
        model_.stabilizeStateValues(stateValues);
    }

    void setIOEvaluationModelTime(const scalar modelTime) const
    {
        model_.setIOEvaluationModelTime(modelTime);
    }

    void debugPrintFields
    (
        const label cellI,
        const scalar t1,
        const scalar t2,
        const scalar step
    ) const
    {
        model_.debugPrintFields(cellI, t1, t2, step);
    }

    void markTransientAlgebraicsDirty() const
    {
        model_.transientAlgebraicsSynchronized_ = false;
    }
};


template<class Hooks, class FieldBuffer>
class BatchedKernelModelPolicy
{
    const Hooks& hooks_;

    template<class CellScratch, class KernelConfig>
    bool rushLarsenEnabledForState
    (
        const label stateI,
        const CellScratch& scratch,
        const KernelConfig& kernel,
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
          ? hooks_.rushLarsenParametersFromHotPathSupport
            (
                stateI,
                scratch.stateValues,
                scratch.rateValues,
                scratch.supportValues,
                steadyState,
                tau
            )
          : hooks_.rushLarsenParameters
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

public:

    explicit BatchedKernelModelPolicy(const Hooks& hooks)
    :
        hooks_(hooks)
    {}

    label nSubsteps() const
    {
        return hooks_.nSubsteps();
    }

    scalar vmInputToState(const scalar vm) const
    {
        return hooks_.vmInputToState(vm);
    }

    template<class CellScratch, class KernelConfig>
    void evaluateScratchAtTime
    (
        const scalar modelTime,
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        if (kernel.useSupportOnlyHotPath)
        {
            hooks_.evaluateHotPathState
            (
                modelTime,
                scratch.stateValues,
                scratch.rateValues,
                scratch.supportValues
            );
        }
        else
        {
            hooks_.evaluateState
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

    template<class CellScratch, class KernelConfig>
    void evaluatePredictorScratchAtTime
    (
        const scalar modelTime,
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        if (kernel.useSupportOnlyHotPath)
        {
            hooks_.evaluateHotPathState
            (
                modelTime,
                scratch.predictorStateValues,
                scratch.predictorRateValues,
                scratch.predictorSupportValues
            );
        }
        else
        {
            hooks_.evaluateState
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

    template<class CellScratch, class KernelConfig>
    void applyExplicitEulerStep
    (
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        for (label stateI = 0; stateI < scratch.stateValues.size(); ++stateI)
        {
            if (!kernel.solveVm && stateI == kernel.vmStateI)
            {
                continue;
            }

            scratch.stateValues[stateI] +=
                kernel.dtSubstep*scratch.rateValues[stateI];
        }
    }

    template<class CellScratch, class KernelConfig>
    void buildPredictorState
    (
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        forAll(scratch.stateValues, stateI)
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

    template<class CellScratch, class KernelConfig>
    void applyCorrectorState
    (
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        forAll(scratch.stateValues, stateI)
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

    template<class CellScratch, class KernelConfig>
    void syncEvaluatedOutputs
    (
        const label cellI,
        const CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        if (kernel.useSupportOnlyHotPath)
        {
            hooks_.scatterCellRates(cellI, scratch.rateValues);
            hooks_.scatterCellSupport(cellI, scratch.supportValues);
        }
        else
        {
            hooks_.scatterCellRatesAndAlgebraics
            (
                cellI,
                scratch.rateValues,
                scratch.algebraicValues
            );
        }
    }

    template<class CellScratch, class KernelConfig>
    scalar ionicCurrentFromScratch
    (
        const CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        if (kernel.useSupportOnlyHotPath)
        {
            return hooks_.ionicCurrentFromHotPathSupport(scratch.supportValues);
        }

        return hooks_.ionicCurrentFromEvaluation(scratch.algebraicValues);
    }

    void stabilizeCellState(const label cellI) const
    {
        hooks_.stabilizeCellState(cellI);
    }

    void stabilizeStateValues(scalarUList& stateValues) const
    {
        hooks_.stabilizeStateValues(stateValues);
    }

    void setIOEvaluationModelTime(const scalar modelTime) const
    {
        hooks_.setIOEvaluationModelTime(modelTime);
    }

    void debugPrintFields
    (
        const label cellI,
        const scalar t1,
        const scalar t2,
        const scalar step
    ) const
    {
        hooks_.debugPrintFields(cellI, t1, t2, step);
    }

    void markTransientAlgebraicsDirty() const
    {
        hooks_.markTransientAlgebraicsDirty();
    }
};


template<class StorageBackend, class ModelPolicy, class FieldBuffer, class StepBuffer>
class BatchedKernelHostBackend
{
    const StorageBackend& storage_;
    const ModelPolicy& policy_;

public:

    BatchedKernelHostBackend
    (
        const StorageBackend& storage,
        const ModelPolicy& policy
    )
    :
        storage_(storage),
        policy_(policy)
    {}

    label nCells() const
    {
        return storage_.nCells();
    }

    label nSubsteps() const
    {
        return policy_.nSubsteps();
    }

    scalar vmInputToState(const scalar vm) const
    {
        return policy_.vmInputToState(vm);
    }

    scalar& stateRef(const label cellI, const label stateI) const
    {
        return storage_.stateRef(cellI, stateI);
    }

    void gatherCellState
    (
        const label cellI,
        FieldBuffer& stateValues
    ) const
    {
        storage_.gatherCellState(cellI, stateValues);
    }

    void scatterCellState
    (
        const label cellI,
        const FieldBuffer& stateValues
    ) const
    {
        storage_.scatterCellState(cellI, stateValues);
    }

    template<class CellScratch, class KernelConfig>
    void evaluateScratchAtTime
    (
        const scalar modelTime,
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        policy_.evaluateScratchAtTime(modelTime, scratch, kernel);
    }

    template<class CellScratch, class KernelConfig>
    void evaluatePredictorScratchAtTime
    (
        const scalar modelTime,
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        policy_.evaluatePredictorScratchAtTime(modelTime, scratch, kernel);
    }

    template<class CellScratch, class KernelConfig>
    void applyExplicitEulerStep
    (
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        policy_.applyExplicitEulerStep(scratch, kernel);
    }

    template<class CellScratch, class KernelConfig>
    void buildPredictorState
    (
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        policy_.buildPredictorState(scratch, kernel);
    }

    template<class CellScratch, class KernelConfig>
    void applyCorrectorState
    (
        CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        policy_.applyCorrectorState(scratch, kernel);
    }

    template<class CellScratch, class KernelConfig>
    void syncEvaluatedOutputs
    (
        const label cellI,
        const CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        policy_.syncEvaluatedOutputs(cellI, scratch, kernel);
    }

    template<class CellScratch, class KernelConfig>
    scalar ionicCurrentFromScratch
    (
        const CellScratch& scratch,
        const KernelConfig& kernel
    ) const
    {
        return policy_.ionicCurrentFromScratch(scratch, kernel);
    }

    void stabilizeCellState(const label cellI) const
    {
        policy_.stabilizeCellState(cellI);
    }

    void stabilizeStateValues(scalarUList& stateValues) const
    {
        policy_.stabilizeStateValues(stateValues);
    }

    StepBuffer& stepField() const
    {
        return storage_.stepField();
    }

    void setIOEvaluationModelTime(const scalar modelTime) const
    {
        policy_.setIOEvaluationModelTime(modelTime);
    }

    void debugPrintFields
    (
        const label cellI,
        const scalar t1,
        const scalar t2,
        const scalar step
    ) const
    {
        policy_.debugPrintFields(cellI, t1, t2, step);
    }

    void markTransientAlgebraicsDirty() const
    {
        policy_.markTransientAlgebraicsDirty();
    }
};


template<class Backend, class KernelConfig, class CellScratch, class SolveProfileTotals>
class BatchedKernelExecutor
{
    using clock = std::chrono::steady_clock;

    const Backend& backend_;
    const KernelConfig& kernel_;
    const bool parallelCellLoops_;
    const bool profileSolve_;
    PtrList<CellScratch>& threadScratch_;
    List<SolveProfileTotals>& threadProfiles_;

    static label cellLoopThreadIndex()
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

public:

    BatchedKernelExecutor
    (
        const Backend& backend,
        const KernelConfig& kernel,
        const bool parallelCellLoops,
        const bool profileSolve,
        PtrList<CellScratch>& threadScratch,
        List<SolveProfileTotals>& threadProfiles
    )
    :
        backend_(backend),
        kernel_(kernel),
        parallelCellLoops_(parallelCellLoops),
        profileSolve_(profileSolve),
        threadScratch_(threadScratch),
        threadProfiles_(threadProfiles)
    {}

    clock::time_point profileStart() const
    {
        return profileSolve_ ? clock::now() : clock::time_point();
    }

    void addProfileTime
    (
        scalar& total,
        const clock::time_point& start
    ) const
    {
        if (profileSolve_)
        {
            total +=
                std::chrono::duration<scalar>(clock::now() - start).count();
        }
    }

    void recordEvaluateCall(label& total) const
    {
        if (profileSolve_)
        {
            ++total;
        }
    }

    template<class VmBuffer>
    void syncVmInputs(const VmBuffer& Vm) const
    {
#ifdef _OPENMP
#pragma omp parallel for if (parallelCellLoops_)
#endif
        for (label cellI = 0; cellI < backend_.nCells(); ++cellI)
        {
            SolveProfileTotals& profile =
                threadProfiles_[cellLoopThreadIndex()];
            const auto vmSyncStart = profileStart();
            if (!kernel_.solveVm)
            {
                backend_.stateRef(cellI, kernel_.vmStateI) =
                    backend_.vmInputToState(Vm[cellI]);
            }
            addProfileTime(profile.vmInputSyncTime, vmSyncStart);
        }
    }

    void emitDebugSample
    (
        const label sampleCell,
        const scalar tStart,
        const scalar tEnd,
        const scalar dtModel
    ) const
    {
        if (sampleCell < 0 || sampleCell >= backend_.nCells())
        {
            return;
        }

        CellScratch& scratch = threadScratch_[0];
        SolveProfileTotals& profile = threadProfiles_[0];
        auto gatherStart = profileStart();
        backend_.gatherCellState(sampleCell, scratch.stateValues);
        addProfileTime(profile.gatherScatterTime, gatherStart);
        scratch.resetPrimary(kernel_.useSupportOnlyHotPath);

        auto evaluateStart = profileStart();
        backend_.evaluateScratchAtTime(tStart, scratch, kernel_);
        addProfileTime(profile.evaluateStateTime, evaluateStart);
        recordEvaluateCall(profile.evaluateStateCalls);

        auto scatterStart = profileStart();
        backend_.syncEvaluatedOutputs(sampleCell, scratch, kernel_);
        addProfileTime(profile.finalSyncTime, scatterStart);
        backend_.stepField()[sampleCell] = dtModel;
        backend_.setIOEvaluationModelTime(tStart);
        backend_.debugPrintFields(sampleCell, tStart, tEnd, dtModel);
        backend_.markTransientAlgebraicsDirty();
    }

    template<class ImBuffer>
    void advanceSubsteps
    (
        const scalar tStart,
        const scalar tEnd,
        ImBuffer& Im
    ) const
    {
        const label nSubsteps = backend_.nSubsteps();
#ifdef _OPENMP
#pragma omp parallel for if (parallelCellLoops_)
#endif
        for (label cellI = 0; cellI < backend_.nCells(); ++cellI)
        {
            CellScratch& scratch = threadScratch_[cellLoopThreadIndex()];
            SolveProfileTotals& profile =
                threadProfiles_[cellLoopThreadIndex()];
            auto gatherStart = profileStart();
            backend_.gatherCellState(cellI, scratch.stateValues);
            addProfileTime(profile.gatherScatterTime, gatherStart);

            for (label substepI = 0; substepI < nSubsteps; ++substepI)
            {
                const scalar substepTime =
                    tStart + scalar(substepI)*kernel_.dtSubstep;

                scratch.resetPrimary(kernel_.useSupportOnlyHotPath);
                scratch.resetPredictor
                (
                    kernel_.useHeun,
                    kernel_.useSupportOnlyHotPath
                );

                auto evaluateStart = profileStart();
                backend_.evaluateScratchAtTime(substepTime, scratch, kernel_);
                addProfileTime(profile.evaluateStateTime, evaluateStart);
                recordEvaluateCall(profile.evaluateStateCalls);
                backend_.stepField()[cellI] = kernel_.dtSubstep;

                auto updateStart = profileStart();
                if (kernel_.useExplicitEuler)
                {
                    backend_.applyExplicitEulerStep(scratch, kernel_);
                    addProfileTime(profile.stateUpdateTime, updateStart);
                }
                else
                {
                    backend_.buildPredictorState(scratch, kernel_);
                    addProfileTime(profile.stateUpdateTime, updateStart);

                    if (kernel_.useHeun)
                    {
                        evaluateStart = profileStart();
                        backend_.evaluatePredictorScratchAtTime
                        (
                            substepTime + kernel_.dtSubstep,
                            scratch,
                            kernel_
                        );
                        addProfileTime(profile.evaluateStateTime, evaluateStart);
                        recordEvaluateCall(profile.evaluateStateCalls);
                    }

                    updateStart = profileStart();
                    backend_.applyCorrectorState(scratch, kernel_);
                    addProfileTime(profile.stateUpdateTime, updateStart);
                }

                auto stabilizeStart = profileStart();
                backend_.stabilizeStateValues(scratch.stateValues);
                addProfileTime(profile.stabilizeTime, stabilizeStart);
            }

            scratch.resetPrimary(kernel_.useSupportOnlyHotPath);
            auto evaluateStart = profileStart();
            backend_.evaluateScratchAtTime(tEnd, scratch, kernel_);
            addProfileTime(profile.evaluateStateTime, evaluateStart);
            recordEvaluateCall(profile.evaluateStateCalls);

            auto finalSyncStart = profileStart();
            backend_.syncEvaluatedOutputs(cellI, scratch, kernel_);
            Im[cellI] = backend_.ionicCurrentFromScratch(scratch, kernel_);
            addProfileTime(profile.finalSyncTime, finalSyncStart);

            auto scatterStart = profileStart();
            backend_.scatterCellState(cellI, scratch.stateValues);
            addProfileTime(profile.gatherScatterTime, scatterStart);
        }
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
