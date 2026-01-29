/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#ifndef activeTensionGenericModel_H
#define activeTensionGenericModel_H

#include "activeTensionModel.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "PtrList.H"
#include "scalarField.H"
#include "dictionary.H"
#include "activeTensionGenericModel_YYYYNames.H"

namespace Foam
{

class activeTensionGenericModel
:
    public activeTensionModel,
    public ODESystem
{
    // Disallow copy
    activeTensionGenericModel(const activeTensionGenericModel&) = delete;
    void operator=(const activeTensionGenericModel&) = delete;

private:

    autoPtr<ODESolver> odeSolver_;
    PtrList<scalarField> internalVariables_;
    PtrList<scalarField> algebraicVariables_;
    scalarField constants_;

public:

    TypeName("activeTensionGenericModel");

    //- Construct from dictionary and number of integration points
    activeTensionGenericModel(const dictionary& dict, const label num);

    virtual ~activeTensionGenericModel() = default;

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
};

} // End namespace Foam

#endif

// ************************************************************************* //
