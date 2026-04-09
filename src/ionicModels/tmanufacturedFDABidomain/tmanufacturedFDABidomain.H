/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef tmanufacturedFDABidomain_H
#define tmanufacturedFDABidomain_H

#include "ionicModel.H"
#include "tmanufacturedFDABidomain_2014Names.H"
#include "OFstream.H"
#include "volFields.H"
#include "overrideTypeName.H"

extern const char* tmanufacturedFDABidomainSTATES_NAMES[BIDOMAIN_NUM_STATES];
extern const char* tmanufacturedFDABidomainALGEBRAIC_NAMES[BIDOMAIN_NUM_ALGEBRAIC];
extern const char* tmanufacturedFDABidomainCONSTANTS_NAMES[BIDOMAIN_NUM_CONSTANTS];

namespace Foam
{

class tmanufacturedFDABidomain
:
    public ionicModel
{
    mutable PtrList<scalarField> STATES_;
    mutable scalarField CONSTANTS_;
    mutable PtrList<scalarField> ALGEBRAIC_;
    mutable PtrList<scalarField> RATES_;
    scalar k_;

    tmanufacturedFDABidomain(const tmanufacturedFDABidomain&) = delete;
    void operator=(const tmanufacturedFDABidomain&) = delete;

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
        return BIDOMAIN_NUM_CONSTANTS;
    }

    virtual const char* const* ioConstantNames() const override
    {
        return tmanufacturedFDABidomainCONSTANTS_NAMES;
    }

    virtual const char* const* ioStateNames() const override
    {
        return tmanufacturedFDABidomainSTATES_NAMES;
    }

    virtual label ioNumStates() const override
    {
        return BIDOMAIN_NUM_STATES;
    }

    virtual const char* const* ioAlgebraicNames() const override
    {
        return tmanufacturedFDABidomainALGEBRAIC_NAMES;
    }

    virtual label ioNumAlgebraic() const override
    {
        return BIDOMAIN_NUM_ALGEBRAIC;
    }

public:

    OverrideTypeName("tmanufacturedFDABidomain");

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
        return "manufacturedFDABidomain";
    }

    virtual label geometricDimension() const override
    {
        return tissue();
    }

    tmanufacturedFDABidomain
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~tmanufacturedFDABidomain();

    virtual void solveODE
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    ) override;

    virtual label nEqns() const override
    {
        return BIDOMAIN_NUM_STATES;
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
