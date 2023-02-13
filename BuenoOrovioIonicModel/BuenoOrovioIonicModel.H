/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::BuenoOrovioIonicModel

Description
    Solves the muscle electro-activation ionic model from A.Bueno-Orovio, E.M.
    Cherry, F.H.Fenton, Minimal model for human ventricular action potentials in
    tissue, J. Theoret. Biol. 253(3) (2008) 544–560.

    This model is defined by three first-order ODEs in terms of three "gating"
    variables.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef BuenoOrovioIonicModel_H
#define BuenoOrovioIonicModel_H

#include "ODESystem.H"
#include "dictionary.H"
#include "volFields.H"
#include "ODESolver.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BuenoOrovioIonicModel Declaration
\*---------------------------------------------------------------------------*/

class BuenoOrovioIonicModel
:
    public ODESystem
{
    // Private Data

        //- Const reference to the dict
        const dictionary& dict_;

        //- Number of integration points, e.g. number of cells
        //  One ODE system will be created for each integration point
        const label nIntegrationPoints_;

        //- Index for monitor point
        //  Gating variables for this point are written to the monitor file
        const label monitorID_;

        //- Monitor file
        OFstream monitorFile_;

        //- ODE solver
        autoPtr<ODESolver> odeSolver_;

        //- List of gating variables for each integration point
        PtrList<scalarField> gatingVariables_;

        //- List of gating variables derivatives for each integration point
        PtrList<scalarField> gatingDerivatives_;

        //- Temporary variable for storing the dimensionless voltage variable
        //  for the integration point  being calculated
        scalar u_;

        //- Model parameters
        const scalar uO_;
        const scalar uU_;
        const scalar thetaV_;
        const scalar thetaW_;
        const scalar thetaVMinus_;
        const scalar thetaO_;
        const scalar tauV1Minus_;
        const scalar tauV2Minus_;
        const scalar tauVPlus_;
        const scalar tauW1Minus_;
        const scalar tauW2Minus_;
        const scalar kWMinus_;
        const scalar uWMinus_;
        const scalar tauWPlus_;
        const scalar tauFi_;
        const scalar tauO1_;
        const scalar tauO2_;
        const scalar tauSo1_;
        const scalar tauSo2_;
        const scalar kSo_;
        const scalar uSo_;
        const scalar tauS1_;
        const scalar tauS2_;
        const scalar kS_;
        const scalar uS_;
        const scalar tauSi_;
        const scalar tauWInfty_;
        const scalar wInftyStar_;

public:

    //- Runtime type information
    TypeName("BuenoOrovio");


    // Constructors

        //- Construct from a dictionary and the number of integration points
        BuenoOrovioIonicModel(const dictionary& dict, const label num);


    //- Destructor
    virtual ~BuenoOrovioIonicModel()
    {}


    // Member Functions

        //- Solve and ODE system for the start time, the time-step and the
        //  voltage. Then place the result in the integration point result list
        virtual void calculateCurrent
        (
            const scalar stepStartTime,
            const scalar timeStep,
            const scalarField& voltage,
            scalarField& result
        );

        //- Number of ODE's to solve
        inline virtual label nEqns() const
        {
            return 3;
        }

        //- Derivatives of each ODE with respect to the primary unknown in that
        //  equation
        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const;

        //- Jacobian for solving the coupled system of ODEs
        virtual void jacobian
        (
            const scalar t,
            const scalarField& y,
            scalarField& dfdt,
            scalarSquareMatrix& dfdy
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
