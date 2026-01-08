/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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
#include "OFstream.H"
#include "dictionary.H"

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

    // ODE solver
    autoPtr<ODESolver> odeSolver_;

    // Internal variables per integration point:
    // y = [Ta, act]
    PtrList<scalarField> internalVariables_;

    // Model parameter(s)
    scalar kTa_;

public:

    TypeName("NashPanfilov");

    //- Construct from dictionary and number of integration points
    NashPanfilov(const dictionary& dict, const label num);

    virtual ~NashPanfilov() = default;

    // --- activeTensionModel API ---

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

    // --- ODESystem API ---

    virtual label nEqns() const override
    {
        return 1;
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


