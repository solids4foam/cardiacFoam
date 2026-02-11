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
    Foam::BuenoOrovio

Description
    Wrapper class for calling from the cellML ionic model for the Bueno Orovio
    2008 reduced ionic model.

SourceFiles
    BuenoOrovio.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef BuenoOrovio_H
#define BuenoOrovio_H

#include "ionicModel.H"
#include "BuenoOrovio_2008Names.H"
#include "OFstream.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Bueno Orovio Declaration
\*---------------------------------------------------------------------------*/

class BuenoOrovio
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
        BuenoOrovio(const BuenoOrovio&) = delete;

        //- No copy assignment
        void operator=(const BuenoOrovio&) = delete;


        public:

        //- Runtime type information
        TypeName("BuenoOrovio");


        //override the only available tissue in this model
        virtual List<word> supportedTissueTypes() const;
        virtual List<word> supportedDimensions() const
        {
            return {}; //No dimensions flag
        }



        // Constructors

            //- Construct given a dictionary, the number of integration points and a
            //  initial time step for the ODE solver
            BuenoOrovio
            (
                const dictionary& dict,
                const label nIntegrationPoints,
                const scalar initialDeltaT,
                const Switch solveVmWithinODESolver = false
            );


        //- Destructor
        virtual ~BuenoOrovio();


        // Member Functions

            virtual void writeHeader(OFstream& output) const;
            //- Write fields to the given file given the current time
            virtual void write(const scalar t, OFstream& output) const;




            virtual void calculateCurrent
            (
                const scalar stepStartTime,
                const scalar deltaT,
                const scalarField& Vm,
                scalarField& Im,
                Field<Field<scalar>>& states
            );


            virtual void solveODE
            (
                const scalar stepStartTime,
                const scalar deltaT,
                const scalarField& Vm,
                scalarField& Im,
                Field<Field<scalar>>& states
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

            // CouplingSignalProvider
            virtual bool hasSignal(const CouplingSignal s) const;
            virtual scalar signal(const label i, const CouplingSignal s) const;


            virtual void updateStatesOld(const Field<Field<scalar>>& states) const;
            virtual void resetStatesToStatesOld(Field<Field<scalar>>& states) const;

                    // Override to provide variable names for export
            virtual Foam::wordList exportedFieldNames() const;
            virtual wordList debugPrintedNames() const;


            // Override to export the variables
            virtual void exportStates
            (
                const Field<Field<scalar>>&,
                PtrList<volScalarField>&
            );

            // Override to print variables in solver
            void debugPrintFields
            (
                label cellI,
                scalar t1,
                scalar t2,
                scalar step
            ) const;

            virtual void sweepCurrent
            (
                const word& currentName,
                scalar Vmin,
                scalar Vmax,
                label nPts,
                const fileName& outputFile
            ) const;


};

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
