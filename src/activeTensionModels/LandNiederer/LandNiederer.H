/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::LandNiederer

Description
    Original seven-state Land et al. (2017) active-tension model.

SourceFiles
    LandNiederer.C
\*---------------------------------------------------------------------------*/

#ifndef LandNiederer_H
#define LandNiederer_H

#include "activeTensionModel.H"
#include "autoPtr.H"
#include "ODESolver.H"
#include "ODESystem.H"
#include "overrideTypeName.H"
#include "PtrList.H"
#include "scalarField.H"
#include "scalarMatrices.H"

#include "LandNiederer_2017Names.H"

namespace Foam
{

class LandNiederer
:
    public activeTensionModel,
    public ODESystem
{
    LandNiederer(const LandNiederer&) = delete;
    void operator=(const LandNiederer&) = delete;

    autoPtr<ODESolver> odeSolver_;

    PtrList<scalarField> STATES_;
    PtrList<scalarField> ALGEBRAIC_;
    PtrList<scalarField> RATES_;

    mutable scalarField CONSTANTS_;
    scalarField prevLambda_;

    mutable scalar currentLambda_;
    mutable scalar currentLambdaRate_;

    void updateDerivedConstants();

protected:

    virtual const PtrList<scalarField>* ioStatesPtr() const override
    {
        return &STATES_;
    }

    virtual const PtrList<scalarField>* ioAlgebraicPtr() const override
    {
        return &ALGEBRAIC_;
    }

    virtual const PtrList<scalarField>* ioRatesPtr() const override
    {
        return &RATES_;
    }

    virtual const scalarField* ioConstantsPtr() const override
    {
        return &CONSTANTS_;
    }

    virtual label ioNumConstants() const override
    {
        return NUM_CONSTANTS;
    }

    virtual const char* const* ioConstantNames() const override;

    virtual const char* const* ioStateNames() const override;

    virtual label ioNumStates() const override
    {
        return NUM_STATES;
    }

    virtual const char* const* ioAlgebraicNames() const override;

    virtual label ioNumAlgebraic() const override
    {
        return NUM_ALGEBRAIC;
    }

    virtual CouplingSignal driveSignal() const override
    {
        return CouplingSignal::CAI;
    }

    virtual scalar driveSignalScaleFactor() const override
    {
        return 1000.0;
    }

    virtual void solveAtPoint
    (
        label i,
        scalar driveVal,
        scalar lambda,
        scalar& Ta
    ) override;

public:

    OverrideTypeName("LandNiederer");

    LandNiederer(const dictionary& dict, const label nIntegrationPoints);

    virtual ~LandNiederer() = default;

    virtual bool supportsRestartState() const override
    {
        return true;
    }

    virtual bool readRestartState(const fvMesh& mesh) override;

    virtual void writeRestartState(const fvMesh& mesh) const override;

    virtual bool restartTension(scalarField& Ta) const override;

    virtual Requirements requirements() const override;

    virtual void preconditionToRestingState(const scalar restingCai) override;

    virtual label nEqns() const override
    {
        return NUM_STATES;
    }

    virtual void derivatives
    (
        const scalar t,
        const scalarField& y,
        scalarField& dydt
    ) const override;

    virtual void jacobian
    (
        const scalar t,
        const scalarField& y,
        scalarField& dfdt,
        scalarSquareMatrix& dfdy
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
