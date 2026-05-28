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

    The class exposes the same compact-support architecture used by the
    other batched models: a cell-major scalar path, a state-major SoA path,
    Rush-Larsen tau/steady-state dispatch, and a thin
    ToRORd_dynClcompactBatched runtime alias that forces compact support.

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

class ToRORd_dynClBatched
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    ToRORd_dynClBatched(const ToRORd_dynClBatched&) = delete;
    void operator=(const ToRORd_dynClBatched&) = delete;

    bool useSoAEvaluator_;
    bool useCompactSupport_;

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

    virtual scalar ionicCurrentFromHotPathSupport
    (
        const scalarUList& supportValues
    ) const override final;

    virtual bool usesCompactSupport() const override final
    {
        return useCompactSupport_;
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
    ) const override final;

    virtual void evaluateHotPathState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& supportValues
    ) const override final;

    virtual bool rushLarsenParametersFromHotPathSupport
    (
        const label stateI,
        const scalarUList& stateValues,
        const scalarUList& rateValues,
        const scalarUList& supportValues,
        scalar& steadyState,
        scalar& tau
    ) const override final;

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

class ToRORd_dynClcompactBatched final
:
    public ToRORd_dynClBatched
{
    ToRORd_dynClcompactBatched(const ToRORd_dynClcompactBatched&) = delete;
    void operator=(const ToRORd_dynClcompactBatched&) = delete;

public:

    OverrideTypeName("ToRORd_dynClcompactBatched");

    ToRORd_dynClcompactBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~ToRORd_dynClcompactBatched() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
