/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Class
    Foam::StewartBatched

Description
    Batched GPU-oriented wrapper for the Stewart 2009 ionic model.

SourceFiles
    StewartBatched.C
\*---------------------------------------------------------------------------*/

#ifndef StewartBatched_H
#define StewartBatched_H

#include "configuredBatchedIonicModel.H"
#include "Stewart_2009Names.H"
#include "StimulusProtocolPOD.H"
#include "overrideTypeName.H"
#include "volFields.H"

#ifdef HAS_CUDA
#include "batchedIonicCoreCUDA.H"
#endif

namespace Foam
{

class StewartBatched final
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    StewartBatched(const StewartBatched&) = delete;
    void operator=(const StewartBatched&) = delete;

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
#endif

protected:

    virtual scalar ionicCurrentFromEvaluation
    (
        const scalarUList& algebraicValues
    ) const override final
    {
        return algebraicValues[Iion_cm];
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

public:

    OverrideTypeName("StewartBatched");

    virtual List<word> supportedTissueTypes() const override final;
    virtual List<word> supportedDimensions() const override final
    {
        return {};
    }

    StewartBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~StewartBatched();

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
