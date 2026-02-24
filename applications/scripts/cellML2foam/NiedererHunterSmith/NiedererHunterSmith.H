/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#ifndef NiedererHunterSmith_H
#define NiedererHunterSmith_H

#include "activeTensionModel.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "PtrList.H"
#include "scalarField.H"
#include "dictionary.H"
#include "NiedererHunterSmith_2006Names.H"

namespace Foam
{

class NiedererHunterSmith
:
    public activeTensionModel,
    public ODESystem
{
    // Disallow copy
    NiedererHunterSmith(const NiedererHunterSmith&) = delete;
    void operator=(const NiedererHunterSmith&) = delete;

private:

    autoPtr<ODESolver> odeSolver_;
    PtrList<scalarField> STATES_;
    PtrList<scalarField> ALGEBRAIC_;
    PtrList<scalarField> RATES_;
    mutable scalarField CONSTANTS_;

public:

    TypeName("NiedererHunterSmith");

    //- Construct from dictionary and number of integration points
    NiedererHunterSmith(const dictionary& dict, const label num);

    virtual ~NiedererHunterSmith() = default;

    virtual Requirements requirements() const override
    {
        Requirements r;
        // Update these flags to match what your model needs (Act, Vm, Ca).
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
