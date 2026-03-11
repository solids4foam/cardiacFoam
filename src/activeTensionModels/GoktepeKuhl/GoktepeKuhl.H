/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#ifndef GoktepeKuhl_H
#define GoktepeKuhl_H

#include "activeTensionModel.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "PtrList.H"
#include "scalarField.H"
#include "dictionary.H"
#include "GoktepeKuhl_2004Names.H"

namespace Foam
{

class GoktepeKuhl
:
    public activeTensionModel,
    public ODESystem
{
    GoktepeKuhl(const GoktepeKuhl&) = delete;
    void operator=(const GoktepeKuhl&) = delete;

private:

    autoPtr<ODESolver> odeSolver_;

    PtrList<scalarField> STATES_;
    PtrList<scalarField> ALGEBRAIC_;
    PtrList<scalarField> RATES_;

    mutable scalarField CONSTANTS_;

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

    //- Returns the configured drive signal (Vm only).
    virtual CouplingSignal driveSignal() const override
    {
        return CouplingSignal::Vm;
    }

    virtual void solveAtPoint
    (
        label i,
        scalar driveVal,
        scalar lambda,
        scalar& Ta
    ) override;

public:

    TypeName("GoktepeKuhl");

    GoktepeKuhl(const dictionary& dict, const label num);

    virtual ~GoktepeKuhl() = default;

    virtual Requirements requirements() const override
    {
        Requirements r;
        r.needVm = true;
        return r;
    }

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
