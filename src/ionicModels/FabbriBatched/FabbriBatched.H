/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Class
    Foam::FabbriBatched

Description
    Batched GPU-oriented wrapper for the Fabbri 2017 ionic model.

SourceFiles
    FabbriBatched.C
\*---------------------------------------------------------------------------*/

#ifndef FabbriBatched_H
#define FabbriBatched_H

#include "configuredBatchedIonicModel.H"
#include "Fabbri_2017Names.H"
#include "StimulusProtocolPOD.H"
#include "overrideTypeName.H"
#include "volFields.H"

#ifdef HAS_CUDA
#include "batchedIonicCoreCUDA.H"
#endif

namespace Foam
{

class FabbriBatched
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    FabbriBatched(const FabbriBatched&) = delete;
    void operator=(const FabbriBatched&) = delete;

protected:

    mutable StimulusProtocolPOD stimulusPOD_;

#ifdef HAS_CUDA
    bool useDevice_;
    mutable batchedIonicCoreCUDA cuda_;
#endif

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
        return true;
    }

    virtual void evaluateState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& algebraicValues
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

    OverrideTypeName("FabbriBatched");

    virtual List<word> supportedTissueTypes() const override final;
    virtual List<word> supportedDimensions() const override final
    {
        return {};
    }

    FabbriBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~FabbriBatched();

    virtual void solveODE
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    ) override;

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

class FabbricompactBatched final
:
    public FabbriBatched
{
    FabbricompactBatched(const FabbricompactBatched&) = delete;
    void operator=(const FabbricompactBatched&) = delete;

public:

    OverrideTypeName("FabbricompactBatched");

    FabbricompactBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~FabbricompactBatched() = default;

    virtual void solveODE
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    ) override;

protected:

#ifdef HAS_CUDA
    void solveOnDevice
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    );
#endif
};

} // End namespace Foam

#endif

// ************************************************************************* //
