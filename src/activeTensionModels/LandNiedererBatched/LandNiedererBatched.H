/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::LandNiedererBatched

Description
    Batched original seven-state Land et al. (2017) active-tension model.
\*---------------------------------------------------------------------------*/

#ifndef LandNiedererBatched_H
#define LandNiedererBatched_H

#include "batchedActiveTensionModel.H"
#include "overrideTypeName.H"

#include "../LandNiederer/LandNiederer_2017Names.H"

namespace Foam
{

class LandNiedererBatched
:
    public batchedActiveTensionModel
{
    LandNiedererBatched(const LandNiedererBatched&) = delete;
    void operator=(const LandNiedererBatched&) = delete;

    mutable scalarField CONSTANTS_;
    mutable scalarField prevLambda_;
    mutable scalarField lambdaRate_;
    mutable scalarField restartTa_;

    void updateDerivedConstants();

protected:

    virtual const char* const* ioConstantNames() const override;
    virtual const char* const* ioStateNames() const override;
    virtual const char* const* ioAlgebraicNames() const override;

    virtual const scalarField* ioConstantsPtr() const override
    {
        return &CONSTANTS_;
    }

    virtual label ioNumConstants() const override
    {
        return NUM_CONSTANTS;
    }

    virtual CouplingSignal driveSignal() const override
    {
        return CouplingSignal::CAI;
    }

    virtual scalar driveSignalScaleFactor() const override
    {
        return 1000.0;
    }

    virtual void evaluateHotPathStateForCell
    (
        label cellI,
        scalar modelTime,
        scalar driveSignal,
        scalar lambda,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& algebraicValues
    ) const override;

    virtual scalar activeTensionFromHotPath
    (
        label cellI,
        const scalarUList& stateValues,
        const scalarUList& algebraicValues
    ) const override
    {
        return algebraicValues[AV_Ta];
    }

public:

    OverrideTypeName("LandNiedererBatched");

    LandNiedererBatched(const dictionary& dict, const label num);

    virtual ~LandNiedererBatched() = default;

    virtual void preconditionToRestingState(const scalar restingCai) override;

    virtual Requirements requirements() const override
    {
        Requirements requirements;
        requirements.needCai = true;
        requirements.needsLambda = true;
        return requirements;
    }

    virtual scalar timeScaleFactor() const override
    {
        return 1000.0;
    }

    virtual void calculateTension
    (
        const scalar t,
        const scalar dt,
        const scalarField& lambda,
        scalarField& Ta
    ) override;

    virtual bool readRestartState(const fvMesh& mesh) override;

    virtual void writeRestartState(const fvMesh& mesh) const override;

    virtual void refreshRestartState(const fvMesh& mesh) override;

    virtual bool restartTension(scalarField& Ta) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
