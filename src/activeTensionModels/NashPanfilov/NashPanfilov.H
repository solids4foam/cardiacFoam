/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#ifndef NashPanfilov_H
#define NashPanfilov_H

#include "activeTensionModel.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "PtrList.H"
#include "scalarField.H"
#include "dictionary.H"
#include "NashPanfilov_2004Names.H"

namespace Foam
{

class NashPanfilov
:
    public activeTensionModel,
    public ODESystem
{
    // Disallow copy
    NashPanfilov(const NashPanfilov&) = delete;
    void operator=(const NashPanfilov&) = delete;

private:

    autoPtr<ODESolver> odeSolver_;

    // States, algebraics, and rates per integration point (ionic-style layout)
    PtrList<scalarField> STATES_;
    PtrList<scalarField> ALGEBRAIC_;
    PtrList<scalarField> RATES_;

    mutable scalarField CONSTANTS_;
    scalar currentAct_;

public:

    TypeName("NashPanfilov");

    //- Construct from dictionary and number of integration points
    NashPanfilov(const dictionary& dict, const label num);

    virtual ~NashPanfilov() = default;

    virtual Requirements requirements() const override
    {
        Requirements r;
        r.needAct = true;
        return r;
    }

    virtual void calculateTension
    (
        const scalar t,
        const scalar dt,
        const scalarField& Act,
        scalarField& Ta
    ) override;

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

    virtual wordList availableFieldNames() const override;
    virtual void fieldValues(const label i, scalarField& values) const override;

    void debugPrintFields
    (
        const label cellI,
        const scalar t1,
        const scalar t2,
        const scalar step,
        const scalar act
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
