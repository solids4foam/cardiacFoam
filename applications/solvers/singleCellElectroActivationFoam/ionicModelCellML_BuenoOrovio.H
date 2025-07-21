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
    Foam::ionicModelCellML

Description
    Wrapper class for calling electrophysiology ionic models from the cellML
    repository.

    Currently the code is setup to run the ten Tusscher, Noble, Noble, panfilov
    (2004) model.

SourceFiles
    ionicModelCellML.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef ionicModelCellML_H
#define ionicModelCellML_H
#include "ionicModelCellML.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ionicModelCellML Declaration
\*---------------------------------------------------------------------------*/

class ionicModelCellML
:
    public ODESystem
{
    // Private Data

        //- ODE solver
        autoPtr<ODESolver> odeSolver_;

        //- States
        mutable scalarField STATES_;

        //- Constants
        mutable scalarField CONSTANTS_;

        //- Algebraic
        mutable scalarField ALGEBRAIC_;

        //- Rates
        mutable scalarField RATES_;

        //- Step size for ODE solver
        mutable scalar step_;

        //flag for tissue type
        mutable scalar tissue_;


    // Private Member Functions

        //- No copy construct
        ionicModelCellML(const ionicModelCellML&) = delete;

        //- No copy assignment
        void operator=(const ionicModelCellML&) = delete;


public:

    //- Runtime type information
    TypeName("ionicModelCellML");

    // Constructors

        //- Construct given a dictionary
        ionicModelCellML(const dictionary& dict);


    //- Destructor
    virtual ~ionicModelCellML();


    // Member Functions

        //- Solve the ODE system from the time = tOld to time = tOld + deltaT
        virtual void solve
        (
            const scalar tOld, const scalar deltaT
        ) const;

        //- Update ALGEBRAIC for the given time
        virtual void computeVariables(const scalar t) const;

        void writeHeader(OFstream& output) const;
        //- Write fields to the given file given the current time
        void write(const scalar t, OFstream& output) const;

        //- Number of ODE's to solve
        inline virtual label nEqns() const
        {
            return 4;
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
