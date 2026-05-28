/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Class
    Foam::PerisYagueBatched

Description
    Batched GPU-oriented wrapper for the Peris-Yague 2022 ionic model.

SourceFiles
    PerisYagueBatched.C
\*---------------------------------------------------------------------------*/

#ifndef PerisYagueBatched_H
#define PerisYagueBatched_H

#include "configuredBatchedIonicModel.H"
#include "PerisYague_2022Names.H"
#include "PerisYague_2022Batch.H"
#include "StimulusProtocolPOD.H"
#include "overrideTypeName.H"
#include "volFields.H"

#ifdef HAS_CUDA
#include "batchedIonicCoreCUDA.H"
#endif

namespace Foam
{

class PerisYagueBatched
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    PerisYagueBatched(const PerisYagueBatched&) = delete;
    void operator=(const PerisYagueBatched&) = delete;

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

    OverrideTypeName("PerisYagueBatched");

    virtual List<word> supportedTissueTypes() const override final;
    virtual List<word> supportedDimensions() const override final
    {
        return {};
    }

    PerisYagueBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~PerisYagueBatched();

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

class PerisYaguecompactBatched final
:
    public PerisYagueBatched
{
    PerisYaguecompactBatched(const PerisYaguecompactBatched&) = delete;
    void operator=(const PerisYaguecompactBatched&) = delete;

public:

    OverrideTypeName("PerisYaguecompactBatched");

    PerisYaguecompactBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~PerisYaguecompactBatched() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
