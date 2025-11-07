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
    Foam::ionicModelFDA

Description
    Abstract base class for ionic models.

SourceFiles
    ionicModelFDA.C

Authors
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef ionicModelFDA_H
#define ionicModelFDA_H

#include "ionicModelFDA.H"
#include "runTimeSelectionTables.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "word.H"
#include "Switch.H"
#include "OFstream.H"
#include "volFields.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class ionicModel Declaration
\*---------------------------------------------------------------------------*/

class ionicModelFDA
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

        //- Flag for tissue type
        mutable label tissue_;

        //- Flag for enable the update of the transmembrane voltage within the
        //  the ionic model ODE solver. For example, if calling for a single
        //  cell, then the governing voltage PDE becomes an ODE so you would
        //  enable this flag
        mutable Switch solveVmWithinODESolver_;


    // Private Member Functions

        //- No copy construct
        ionicModelFDA(const ionicModelFDA&) = delete;

        //- No copy assignment
        void operator=(const ionicModelFDA&) = delete;


        



        
public:

    void setTissueFromDict();

    //- Runtime type information
    TypeName("ionicModelFDA");


    // Declare run-time constructor selection table

        //- Ionic model run-time selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            ionicModelFDA,
            dictionary,
            (
                const dictionary& dict,
                const label nIntegrationPoints,
                const scalar initialDeltaT,
                const Switch solveVmWithinODESolver
            ),
            (dict, nIntegrationPoints, initialDeltaT, solveVmWithinODESolver)
        );

    // Selectors

        //- Create a ionic model
        static autoPtr<ionicModelFDA> New
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );

    // Constructors

        //- Construct given a dictionary, the number of integration points and a
        //  initial time step for the ODE solver
        ionicModelFDA
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );


    //- Destructor
    virtual ~ionicModelFDA();


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

        //- Getter for tissue flag
        const Switch& solveVmWithinODESolver() const
        {
            return solveVmWithinODESolver_;
        }

        //- Solve and ODE system for the start time, the time-step and the
        //  voltage. Then place the result in the integration point result list

        virtual void calculateCurrent
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            scalarField& u1m,
            scalarField& u2m,
            scalarField& u3m
        )
        {
            notImplemented("WIP");
        }
        
        virtual void calculateGating
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            scalarField& u1m,
            scalarField& u2m,
            scalarField& u3m
        )
        {
            notImplemented("WIP");
        }

        virtual void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            scalarField& u1m,
            scalarField& u2m,
            scalarField& u3m
        )
        {
            notImplemented("WIP");
        }
        
        virtual void initializeFields
        (
            volScalarField& Vm,
            volScalarField& u1,
            volScalarField& u2,
            volScalarField& u3,
            const volVectorField& C
        )
        {
            notImplemented("WIP");
        }


        //Supported words by the solver 
        virtual List<word> supportedTissues() const
        {
            return {"1D", "2D", "3D"};
        }

        virtual void writeHeader(OFstream& os) const 
        {

        }
        virtual void write(const scalar t, OFstream& os) const 
        {

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
    };  

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
