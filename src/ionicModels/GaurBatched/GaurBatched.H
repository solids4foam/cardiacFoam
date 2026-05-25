/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::GaurBatched

Description
    Batched / GPU-portable wrapper for the Gaur 2021 ionic model.

    --- Stage in the GPU port ---

    This is the *first stage* of the Gaur GPU port. It plugs Gaur into
    the existing cell-major batched executor (RL+Heun on host, OpenMP-
    parallel cell loop optional). It is the equivalent of what
    `TNNPBatched` looks like today: a runtime-selectable wrapper around the
    full `GaurcomputeVariables` evaluator, exercised through the
    validated cell-major path.

    A future stage will add a state-major SoA evaluator
    (`Gaur_2021Batch.H`), opt-in `useSoAEvaluator` switch,
    device buffers, and `GaurBatched_cuda.cu` — mirroring exactly what
    `BuenoOrovioBatched` looks like today. That stage is not in this
    commit; see `gpu_validation_strategy.md` for the recipe.

    Until that stage lands, `GaurBatched` is functionally similar to the
    legacy `Gaur` class but uses the host-batched RL+Heun integrator
    instead of OpenFOAM's RKF45. The default integrator
    `batchedIntegrator rlGatesHeunRest;` keeps gate states on Rush-
    Larsen and the rest on Heun — the same accuracy class as the
    legacy CPU model in practice.

SourceFiles
    GaurBatched.C
\*---------------------------------------------------------------------------*/

#ifndef GaurGPU_H
#define GaurGPU_H

#include "configuredBatchedIonicModel.H"
#include "OFstream.H"
#include "volFields.H"
#include "Gaur_2021Names.H"
#include "Gaur_2021Batch.H"
#include "StimulusProtocolPOD.H"
#include "overrideTypeName.H"

#ifdef HAS_CUDA
#include "batchedIonicCoreCUDA.H"
#endif

namespace Foam
{

class GaurBatched final
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    GaurBatched(const GaurBatched&) = delete;
    void operator=(const GaurBatched&) = delete;

    bool useSoAEvaluator_;
    mutable StimulusProtocolPOD stimulusPOD_;

    void solveBatched
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    );

#ifdef HAS_CUDA
    bool useDevice_;
    mutable batchedIonicCoreCUDA cuda_;
#endif

protected:

    virtual scalar ionicCurrentFromEvaluation
    (
        const scalarUList& algebraicValues
    ) const final
    {
        return algebraicValues[Iion_cm];
    }

    virtual void evaluateState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& algebraicValues
    ) const final;

    virtual bool rushLarsenParameters
    (
        const label stateI,
        const scalarUList& stateValues,
        const scalarUList& rateValues,
        const scalarUList& algebraicValues,
        scalar& steadyState,
        scalar& tau
    ) const final;

public:

    OverrideTypeName("GaurBatched");

    virtual List<word> supportedTissueTypes() const override final;
    virtual List<word> supportedDimensions() const override final
    {
        return {};
    }

    GaurBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~GaurBatched();

    virtual void solveODE
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    ) override final;

    virtual void prepareIOAccess
    (
        const wordList& requestedNames,
        const bool needsAlgebraics
    ) const override final;

    virtual void importFields
    (
        const volScalarField& Vm,
        const wordList& fieldNames,
        const PtrList<volScalarField>& inFields
    ) override final;

    virtual void derivatives
    (
        const scalar t,
        const scalarField& y,
        scalarField& dydt
    ) const override final;


};

} // End namespace Foam

#endif

// ************************************************************************* //
