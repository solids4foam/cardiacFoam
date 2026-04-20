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
    Foam::PerisYague_2022

Description
    Wrapper class for the Peris-Yaguë atrial ionic model 2022.

SourceFiles
    PerisYague_2022.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef PerisYague_2022_H
#define PerisYague_2022_H

#include "ionicModel.H"
#include "OFstream.H"
#include "volFields.H"
#include "PerisYague_2022Names.H"
#include "dictionary.H"
#include "overrideTypeName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PerisYague_2022 Declaration
\*---------------------------------------------------------------------------*/

class PerisYague_2022
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
        PerisYague_2022(const PerisYague_2022&) = delete;

        //- No copy assignment
        void operator=(const PerisYague_2022&) = delete;

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
        OverrideTypeName("PerisYague_2022");


        //override the only available tissue in this model
        virtual List<word> supportedTissueTypes() const override;
        virtual List<word> supportedDimensions() const override
        {
            return {}; //No dimensions flag
        }



        // Constructors

            //- Construct given a dictionary, the number of integration points and a
            //  initial time step for the ODE solver
            PerisYague_2022
            (
                const dictionary& dict,
                const label nIntegrationPoints,
                const scalar initialDeltaT,
                const Switch solveVmWithinODESolver = false
            );


        //- Destructor
        virtual ~PerisYague_2022();


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
