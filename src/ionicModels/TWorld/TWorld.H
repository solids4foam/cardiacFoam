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
    Foam::TWorld

Description
    Wrapper class for calling from the cellML ionic model for the TWorld model.

SourceFiles
    TWorld.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef TWorld_H
#define TWorld_H

#include "ionicModel.H"
#include "OFstream.H"
#include "volFields.H"
#include "TWorld_2024Names.H"
#include "overrideTypeName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class TWorld Declaration
\*---------------------------------------------------------------------------*/

class TWorld
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

        //- Per-integration-point constants when transmural heterogeneity is on
        mutable PtrList<scalarField> HETEROGENEOUS_CONSTANTS_;

        //- Active integration point for derivatives() (no index argument)
        mutable label activeIntegrationPoint_ = 0;


    // Private Member Functions

        //- No copy construct
        TWorld(const TWorld&) = delete;

        //- No copy assignment
        void operator=(const TWorld&) = delete;

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

        virtual scalarField* ioMutableConstantsPtr() const override
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
    OverrideTypeName("TWorld");


    //- Tissue types supported by this model (endo / M / epi)
    virtual List<word> supportedTissueTypes() const override;

    //- Sex types supported by this model (neutral / male / female)
    virtual List<word> supportedSexTypes() const override;

    //- Per-point constants accessor (heterogeneous when configured)
    scalarField& constants(const label integrationPtI) const;

    //- Full constants array for a tissue flag (endo 3, M 2, epi 1)
    virtual scalarField constantsForTissue
    (
        const label tissueFlag
    ) const override;

    //- Build per-point transmural heterogeneity from electroProperties
    virtual void configureIonicHeterogeneity
    (
        const scalarField& transmuralDistance,
        const dictionary& heterogeneityDict
    ) override;
    virtual List<word> supportedDimensions() const override
    {
        return {}; //No dimensions flag
    }



    // Constructors

        //- Construct given a dictionary, the number of integration points and a
        TWorld
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );


    //- Destructor
    virtual ~TWorld();


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
