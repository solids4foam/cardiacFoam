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
    Foam::Fabbri

Description
    Wrapper class for calling from the cellML ionic model for the Fabbri model.

SourceFiles
    Fabbri.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef Fabbri_H
#define Fabbri_H

#include "ionicModel.H"
#include "OFstream.H"
#include "volFields.H"
#include "Fabbri_2017Names.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Fabbri Declaration
\*---------------------------------------------------------------------------*/

class Fabbri
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
        Fabbri(const Fabbri&) = delete;

        //- No copy assignment
        void operator=(const Fabbri&) = delete;

protected:
        virtual const PtrList<scalarField>* ioStatesPtr() const override
        {
            return &STATES_;
        }

        virtual const PtrList<scalarField>* ioAlgebraicPtr() const override
        {
            return &ALGEBRAIC_;
        }

        virtual const PtrList<scalarField>* ioRatesPtr() const override
        {
            return &RATES_;
        }

        virtual const scalarField* ioConstantsPtr() const override
        {
            return &CONSTANTS_;
        }

        virtual label ioNumConstants() const override
        {
            return NUM_CONSTANTS;
        }

        virtual const char* const* ioConstantNames() const override;

        virtual const char* const* ioStateNames() const override;

        virtual label ioNumStates() const override
        {
            return NUM_STATES;
        }

        virtual const char* const* ioAlgebraicNames() const override;

        virtual label ioNumAlgebraic() const override
        {
            return NUM_ALGEBRAIC;
        }


        public:

        //- Runtime type information
        TypeName("Fabbri");


        //override the only available tissue in this model
        virtual List<word> supportedTissueTypes() const override;
        virtual List<word> supportedDimensions() const override
        {
            return {}; //No dimensions flag
        }



        // Constructors

            //- Construct given a dictionary, the number of integration points and a
            //  initial time step for the ODE solver
            Fabbri
            (
                const dictionary& dict,
                const label nIntegrationPoints,
                const scalar initialDeltaT,
                const Switch solveVmWithinODESolver = false
            );


        //- Destructor
        virtual ~Fabbri();


        // Member Functions



            virtual void solveODE
            (
                const scalar stepStartTime,
                const scalar deltaT,
                const scalarField& Vm,
                scalarField& Im
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

            virtual void sweepCurrent
            (
                const word& currentName,
                scalar Vmin,
                scalar Vmax,
                label nPts,
                const fileName& outputFile
            ) const override;

            virtual wordList availableSweepCurrents() const override;


};

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
