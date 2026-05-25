/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::ToRORd_dynClBatched

Description
    Batched / GPU-portable wrapper for the ToRORd_dynCl 2023 ionic model
    (Tomek-O'Hara-Rudy with dynamic chloride).

    --- Stage in the GPU port ---

    First stage: cell-major batched executor wrapper. Same architectural
    level as GaurBatched / TNNPBatched today — the model goes through
    the validated cell-major executor (Heun by default for ToRORd, since
    no Rush-Larsen dispatch table is provided here). A future iteration
    may add:

      1. A `rushLarsenParameters` override final exposing the model's many
         (steady-state, tau) gate pairs to the executor's RL+Heun path.
      2. A state-major SoA Batch evaluator (`ToRORd_dynCl_2023Batch.H`)
         + `useSoAEvaluator` switch + `*_cuda.cu` placeholder, mirroring
         BuenoOrovioBatched.

    For now this class exists primarily as a stress-test target for
    measurement work: ~40 states is the largest active model in the
    tree and exercises the integrator + per-cell scratch + write-seam
    machinery much harder than the smaller models.

SourceFiles
    ToRORd_dynClBatched.C
\*---------------------------------------------------------------------------*/

#ifndef ToRORd_dynClBatched_H
#define ToRORd_dynClBatched_H

#include "configuredBatchedIonicModel.H"
#include "ToRORd_dynCl_2023Batch.H"
#include "StimulusProtocolPOD.H"
#include "OFstream.H"
#include "volFields.H"
#include "ToRORd_dynCl_2023Names.H"
#include "overrideTypeName.H"

#ifdef HAS_CUDA
#include "batchedIonicCoreCUDA.H"
#endif

namespace Foam
{

class ToRORd_dynClBatched final
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    ToRORd_dynClBatched(const ToRORd_dynClBatched&) = delete;
    void operator=(const ToRORd_dynClBatched&) = delete;

    bool useSoAEvaluator_;

#ifdef HAS_CUDA
    bool useDevice_;
#endif

    mutable StimulusProtocolPOD stimulusPOD_;

    void solveBatched
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    );

#ifdef HAS_CUDA
    mutable batchedIonicCoreCUDA cuda_;
#endif // HAS_CUDA

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

public:

    OverrideTypeName("ToRORd_dynClBatched");

    virtual List<word> supportedTissueTypes() const override final;
    virtual List<word> supportedDimensions() const override final
    {
        return {};
    }

    ToRORd_dynClBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~ToRORd_dynClBatched();

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
