/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Class
    Foam::CourtemancheBatched

Description
    Batched GPU-oriented wrapper for the Courtemanche 1998 ionic model.

SourceFiles
    CourtemancheBatched.C
\*---------------------------------------------------------------------------*/

#ifndef CourtemancheBatched_H
#define CourtemancheBatched_H

#include "configuredBatchedIonicModel.H"
#include "Courtemanche_1998Names.H"
#include "Courtemanche_1998Batch.H"
#include "StimulusProtocolPOD.H"
#include "overrideTypeName.H"
#include "volFields.H"

#ifdef HAS_CUDA
#include "batchedIonicCoreCUDA.H"
#endif

namespace Foam
{

class CourtemancheBatched
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    CourtemancheBatched(const CourtemancheBatched&) = delete;
    void operator=(const CourtemancheBatched&) = delete;

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
#endif

protected:

    virtual scalar ionicCurrentFromEvaluation
    (
        const scalarUList& algebraicValues
    ) const override final
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
    ) const override final;

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

    OverrideTypeName("CourtemancheBatched");

    virtual List<word> supportedTissueTypes() const override final;
    virtual List<word> supportedDimensions() const override final
    {
        return {};
    }

    CourtemancheBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~CourtemancheBatched();

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

class CourtemanchecompactBatched final
:
    public CourtemancheBatched
{
    CourtemanchecompactBatched(const CourtemanchecompactBatched&) = delete;
    void operator=(const CourtemanchecompactBatched&) = delete;

public:

    OverrideTypeName("CourtemanchecompactBatched");

    CourtemanchecompactBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~CourtemanchecompactBatched() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
