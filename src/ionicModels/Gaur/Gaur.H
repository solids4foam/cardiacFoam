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
    Wrapper class for calling from the cellML ionic model for the Gaur model 2021.

SourceFiles
    Gaur.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef Gaur_H
#define Gaur_H

#include "ionicModel.H"
#include "OFstream.H"
#include "volFields.H"
#include "Gaur_2021Names.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class manufacturedFDA Declaration
\*---------------------------------------------------------------------------*/

class Gaur
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
        Gaur(const Gaur&) = delete;

        //- No copy assignment
        void operator=(const Gaur&) = delete;


public:

    //- Runtime type information
    TypeName("Gaur");

    virtual List<word> supportedTissueTypes() const override;
    virtual List<word> supportedDimensions() const override
    {
        return {}; //No dimensions flag
    }


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

        virtual void writeHeader(OFstream& output) const override;
        //- Write fields to the given file given the current time
        virtual void write(const scalar t, OFstream& output) const override;


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


        virtual void updateStatesOld(const Field<Field<scalar>>& states) const;
        virtual void resetStatesToStatesOld(Field<Field<scalar>>& states) const;

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

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

