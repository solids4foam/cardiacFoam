/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
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
    Foam::ionicModel

Description
    Abstract base class for ionic models.

SourceFiles
    ionicModel.C

Authors
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef ionicModel_H
#define ionicModel_H

#include "ionicModel.H"
#include "runTimeSelectionTables.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "word.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class ionicModel Declaration
\*---------------------------------------------------------------------------*/

class ionicModel
:
    public ODESystem
{
    // Private Data

    // Private Data

        //- ODE solver
        autoPtr<ODESolver> odeSolver_;

        //- Copy of the model dictionary
        dictionary dict_;

        //- Step size for ODE solver
        mutable scalarList step_;

        //flag for tissue type
        mutable label tissue_;

    // Private Member Functions

        //- No copy construct
        ionicModel(const ionicModel&) = delete;

        //- No copy assignment
        void operator=(const ionicModel&) = delete;


public:

    //- Runtime type information
    TypeName("ionicModel");


    // Declare run-time constructor selection table

        //- Ionic model run-time selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            ionicModel,
            dictionary,
            (
                const dictionary& dict,
                const label nIntegrationPoints,
                const scalar initialDeltaT
            ),
            (dict, nIntegrationPoints, initialDeltaT)
        );

    // Selectors

        //- Create a ionic model
        static autoPtr<ionicModel> New
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT
        );

    // Constructors

        //- Construct given a dictionary, the number of integration points and a
        //  initial time step for the ODE solver
        ionicModel
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT
        );


    //- Destructor
    virtual ~ionicModel();


    // Member Functions

        //- ODE solver
        ODESolver& odeSolver()
        {
            if (!odeSolver_)
            {
                odeSolver_.reset(ODESolver::New(*this, dict_));
            }

            return odeSolver_();
        }

        //- Step size for ODE solver
        scalarList& step()
        {
            return step_;
        }

        //- Getter for tissue flag
        label tissue() const
        {
            return tissue_;
        }

        //- Solve and ODE system for the start time, the time-step and the
        //  voltage. Then place the result in the integration point result list
        virtual void calculateCurrent
        (
            const scalar stepStartTime,
            const scalar timeStep,
            const scalarField& voltage,
            scalarField& result
        )
        {
            notImplemented("WIP");
        }

        //- Number of ODE's to solve
        virtual label nEqns() const
        {
            notImplemented("WIP");

            return 1;
        }

        //- Derivatives of each ODE with respect to the primary unknown in that
        //  equation
        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const
        {
            notImplemented("WIP");
        }

        //- Jacobian for solving the coupled system of ODEs
        virtual void jacobian
        (
            const scalar t,
            const scalarField& y,
            scalarField& dfdt,
            scalarSquareMatrix& dfdy
        ) const
        {
            notImplemented
            (
                "Please define this function in the derived class!"
            );
        };

        //This functions are used for the single cell electro-activation solver

        // //- Solve the ODE system from the time = tOld to time = tOld + deltaT
        // virtual void writeHeader(Ostream& os) const = 0;
        // virtual void solve(const scalar tOld, const scalar deltaT) = 0;
        // virtual void computeVariables(const scalar VOI) = 0;
        // virtual void write(const scalar VOI, Ostream& os) const = 0;
        // // inline virtual label nEqns() const
        // {
        //     return NUM_STATES;
        // }
    };  

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
