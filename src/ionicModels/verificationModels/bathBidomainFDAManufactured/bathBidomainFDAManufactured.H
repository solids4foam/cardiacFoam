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
    Foam::bathBidomainFDAManufactured

Description
    Bath-bidomain variant of the FDA manufactured ionic model.

    The bath variant differs from bidomainFDAManufactured by exactly one
    invariant: the analytical state u3 = -alpha*x/se, which carries the
    extracellular surface current along the bath boundary. This invariant
    is baked into the dedicated header bathBidomainFDAManufactured_2014.H;
    the initial volScalarField value of u3 is set by the matching verifier
    (manufacturedFDABathBidomainVerifier::preProcess) and then preserved
    by the ODE which sets RATES[u3] = 0.

    No `bath` runtime switch exists on this model -- the user selects it
    by name. This replaces the legacy `manufacturedBidomain { bath yes; }`
    fork inside bidomainFDAManufactured.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef bathBidomainFDAManufactured_H
#define bathBidomainFDAManufactured_H

#include "ionicModel.H"
#include "bathBidomainFDAManufactured_2014Names.H"
#include "OFstream.H"
#include "volFields.H"
#include "overrideTypeName.H"

extern const char*
    bathBidomainFDAManufacturedSTATES_NAMES[BATH_BIDOMAIN_NUM_STATES];
extern const char*
    bathBidomainFDAManufacturedALGEBRAIC_NAMES[BATH_BIDOMAIN_NUM_ALGEBRAIC];
extern const char*
    bathBidomainFDAManufacturedCONSTANTS_NAMES[BATH_BIDOMAIN_NUM_CONSTANTS];

namespace Foam
{

class bathBidomainFDAManufactured
:
    public ionicModel
{
    mutable PtrList<scalarField> STATES_;
    mutable scalarField CONSTANTS_;
    mutable PtrList<scalarField> ALGEBRAIC_;
    mutable PtrList<scalarField> RATES_;
    scalar k_;

    void refreshAlgebraics() const;

    bathBidomainFDAManufactured(const bathBidomainFDAManufactured&) = delete;
    void operator=(const bathBidomainFDAManufactured&) = delete;

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
        return BATH_BIDOMAIN_NUM_CONSTANTS;
    }

    virtual const char* const* ioConstantNames() const override
    {
        return bathBidomainFDAManufacturedCONSTANTS_NAMES;
    }

    virtual const char* const* ioStateNames() const override
    {
        return bathBidomainFDAManufacturedSTATES_NAMES;
    }

    virtual label ioNumStates() const override
    {
        return BATH_BIDOMAIN_NUM_STATES;
    }

    virtual const char* const* ioAlgebraicNames() const override
    {
        return bathBidomainFDAManufacturedALGEBRAIC_NAMES;
    }

    virtual label ioNumAlgebraic() const override
    {
        return BATH_BIDOMAIN_NUM_ALGEBRAIC;
    }

public:

    OverrideTypeName("bathBidomainFDAManufactured");

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
        return "bathBidomainFDAManufactured";
    }

    virtual label geometricDimension() const override
    {
        return tissue();
    }

    bathBidomainFDAManufactured
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~bathBidomainFDAManufactured();

    virtual void solveODE
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    ) override;

    //- Vm is advanced linearly across the step from the supplied rate
    virtual bool supportsVmRateCoupling() const override
    {
        return true;
    }

    //- Iion can be recomputed at any Vm from the stored state
    virtual bool supportsIonicCurrentEvaluation() const override
    {
        return true;
    }

    virtual void evaluateIonicCurrent
    (
        const scalar t,
        const scalarField& Vm,
        scalarField& Im
    ) override;

    virtual void exportStates
    (
        PtrList<volScalarField>& outFields
    ) override;

    virtual label nEqns() const override
    {
        return BATH_BIDOMAIN_NUM_STATES;
    }

    virtual void derivatives
    (
        const scalar t,
        const scalarField& y,
        scalarField& dydt
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
