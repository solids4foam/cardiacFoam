/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

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

#include "ionicModel.H"
#include "tmanufacturedFDA_2014Names.H"
#include "OFstream.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class manufacturedFDA Declaration
\*---------------------------------------------------------------------------*/

class tmanufacturedFDA
:
    public ionicModel
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
        tmanufacturedFDA(const tmanufacturedFDA&) = delete;

        //- No copy assignment
        void operator=(const tmanufacturedFDA&) = delete;


public:

    //- Runtime type information
    TypeName("tmanufacturedFDA");


    //override the only available tissue in this model
    virtual List<word> supportedDimensions() const override;


    // Constructors

        //- Construct given a dictionary, the number of integration points and a
        //  initial time step for the ODE solver
        tmanufacturedFDA
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );


    //- Destructor
    virtual ~tmanufacturedFDA();


    // Member Functions

        //virtual void writeHeader(OFstream& output) const override;
        //- Write fields to the given file given the current time
        //virtual void write(const scalar t, OFstream& output) const override;

        virtual bool hasManufacturedSolution() const override
        {
            return true;
        }

        virtual void calculateCurrent
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            Field<Field<scalar>>& states
        ) override;

        virtual void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            Field<Field<scalar>>& states
        ) override;

        void initializeFields
        (
            volScalarField& Vm,
            volScalarField& u1,
            volScalarField& u2,
            volScalarField& u3,
            const volVectorField& C
        ) override;


        //- Number of ODE's to solve
        virtual label nEqns() const override
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
        ) const override;

        void updateStatesOld() const;
        void resetStatesToStatesOld() const;
        // Override to provide variable names for export
        virtual Foam::wordList exportedFieldNames() const override;
        virtual wordList debugPrintedNames() const override;


        // Override to export the variables
        virtual void exportStates
        (
            const Field<Field<scalar>>&,
            List<volScalarField*>&
        ) override;

        // Override to print variables in solver
        void debugPrintFields
        (
            label cellI,
            scalar t1,
            scalar t2,
            scalar step
        ) const override;




};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
