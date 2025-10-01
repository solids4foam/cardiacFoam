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
    Foam::Gaur

Description
    Wrapper class for calling from the cellML ionic model for the ten Tusscher,
    Noble, Noble, Panfilov, 2004, model.

SourceFiles
    Gaur.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef Gaur_H
#define Gaur_H

#include "ionicModel.H"
#include "Gaur_2021Names.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Gaur Declaration
\*---------------------------------------------------------------------------*/

class Gaur
:
    public ionicModel
{
    // Private Data

        //- States
        mutable PtrList<scalarField> STATES_;

        //- Constants
        mutable scalarField CONSTANTS_;

        //- Algebraic
        mutable PtrList<scalarField> ALGEBRAIC_;

        //- Rates
        mutable PtrList<scalarField> RATES_;


    // Private Member Functions

        //- No copy construct
        Gaur(const Gaur&) = delete;

        //- No copy assignment
        void operator=(const Gaur&) = delete;


public:

    //- Runtime type information
    TypeName("Gaur");

    // Constructors

        //- Construct given a dictionary, the number of integration points and a
        //  initial time step for the ODE solver
        Gaur
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );


    //- Destructor
    virtual ~Gaur();


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
        virtual label nEqns() const
        {
            return NUM_STATES;
        }

        //- Derivatives of each ODE with respect to the primary unknown in that
        //  equation
        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
