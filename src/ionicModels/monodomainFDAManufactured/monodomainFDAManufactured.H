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
    Foam::monodomainFDAManufactured

Description
    Wrapper class for calling from the cellML ionic model for the ten Tusscher,
    Noble, Noble, Panfilov, 2004, model.

SourceFiles
    monodomainFDAManufactured.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef monodomainFDAManufactured_H
#define monodomainFDAManufactured_H

#include "ionicModel.H"
#include "monodomainFDAManufactured_2014Names.H"
#include "OFstream.H"
#include "volFields.H"
#include "overrideTypeName.H"

extern const char* monodomainFDAManufacturedSTATES_NAMES[NUM_STATES];
extern const char* monodomainFDAManufacturedALGEBRAIC_NAMES[NUM_ALGEBRAIC];
extern const char* monodomainFDAManufacturedCONSTANTS_NAMES[NUM_CONSTANTS];

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class monodomainFDAManufactured Declaration
\*---------------------------------------------------------------------------*/

class monodomainFDAManufactured
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
    monodomainFDAManufactured(const monodomainFDAManufactured&) = delete;

    //- No copy assignment
    void operator=(const monodomainFDAManufactured&) = delete;


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

    virtual const char* const* ioConstantNames() const override
    {
        return monodomainFDAManufacturedCONSTANTS_NAMES;
    }

    virtual const char* const* ioStateNames() const override
    {
        return monodomainFDAManufacturedSTATES_NAMES;
    }

    virtual label ioNumStates() const override
    {
        return NUM_STATES;
    }

    virtual const char* const* ioAlgebraicNames() const override
    {
        return monodomainFDAManufacturedALGEBRAIC_NAMES;
    }

    virtual label ioNumAlgebraic() const override
    {
        return NUM_ALGEBRAIC;
    }


public:

    //- Runtime type information
    OverrideTypeName("monodomainFDAManufactured");

    // Override the only available dimension family in this model
    virtual List<word> supportedDimensions() const override;

    virtual word tissueName() const override
    {
        return (tissue() == 1) ? "1D"
             : (tissue() == 2) ? "2D"
             : (tissue() == 3) ? "3D"
             :                   "unknown";
    }

    virtual word verificationFamily() const override
    {
        return "monodomainFDAManufactured";
    }

    virtual label geometricDimension() const override
    {
        return tissue();
    }


    // Constructors

    //- Construct given a dictionary, the number of integration points and a
    //  initial time step for the ODE solver
    monodomainFDAManufactured
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    //- Destructor
    virtual ~monodomainFDAManufactured();


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
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
