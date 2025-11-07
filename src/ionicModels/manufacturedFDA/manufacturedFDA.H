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
    Foam::manufacturedFDA

Description
    Wrapper class for calling from the cellML ionic model for the ten Tusscher,
    Noble, Noble, Panfilov, 2004, model.

SourceFiles
    manufacturedFDA.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef manufacturedFDA_H
#define manufacturedFDA_H

#include "ionicModelFDA.H"
#include "manufacturedFDA_2014Names.H"
#include "OFstream.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class manufacturedFDA Declaration
\*---------------------------------------------------------------------------*/

class manufacturedFDA
:
    public ionicModelFDA
{
    // Private Data

        //- States
        mutable PtrList<scalarField> STATES_;

        //- States - old time
        mutable PtrList<scalarField> STATES_OLD_;

        //- Constants
        mutable scalarField CONSTANTS_;

        //- Algebraic
        mutable PtrList<scalarField> ALGEBRAIC_;

        //- Rates
        mutable PtrList<scalarField> RATES_;


    // Private Member Functions

        //- No copy construct
        manufacturedFDA(const manufacturedFDA&) = delete;

        //- No copy assignment
        void operator=(const manufacturedFDA&) = delete;


public:

    //- Runtime type information
    TypeName("manufacturedFDA");


    //override the only available tissue in this model
    virtual List<word> supportedTissues() const;

    // Constructors

        //- Construct given a dictionary, the number of integration points and a
        //  initial time step for the ODE solver
        manufacturedFDA
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );


    //- Destructor
    virtual ~manufacturedFDA();


    // Member Functions

        virtual void writeHeader(OFstream& output) const;
        //- Write fields to the given file given the current time
        virtual void write(const scalar t, OFstream& output) const;

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
        );
        virtual void calculateGating
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            scalarField& u1m,
            scalarField& u2m,
            scalarField& u3m
        );
        virtual void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            scalarField& u1m,
            scalarField& u2m,
            scalarField& u3m
        );


        virtual void initializeFields
        (
            volScalarField& Vm,
            volScalarField& u1,
            volScalarField& u2,
            volScalarField& u3,
            const volVectorField& C
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

        //- Reset STATES to STATES_OLD
        virtual void resetStatesToStatesOld() const
        {
            STATES_ = STATES_OLD_;
        }

        //- Update STATES_OLD from STATES
        virtual void updateStatesOld() const
        {
            STATES_OLD_ = STATES_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
