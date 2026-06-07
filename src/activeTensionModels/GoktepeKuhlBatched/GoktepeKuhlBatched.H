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
    Foam::GoktepeKuhlBatched

Description
    Goktepe-Kuhl (2010) active tension model implemented for batched GPU
    execution.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef GoktepeKuhlBatched_H
#define GoktepeKuhlBatched_H

#include "batchedActiveTensionModel.H"
#include "../GoktepeKuhl/GoktepeKuhl_2004Names.H"
#include "overrideTypeName.H"

namespace Foam
{

class GoktepeKuhlBatched
:
    public batchedActiveTensionModel
{
    GoktepeKuhlBatched(const GoktepeKuhlBatched&) = delete;
    void operator=(const GoktepeKuhlBatched&) = delete;

    mutable scalarField CONSTANTS_;

protected:

    virtual const char* const* ioConstantNames() const override;
    virtual const char* const* ioStateNames() const override;
    virtual const char* const* ioAlgebraicNames() const override;

    virtual const scalarField* ioConstantsPtr() const override
    {
        return &CONSTANTS_;
    }

    virtual label ioNumConstants() const override
    {
        return NUM_CONSTANTS;
    }

    virtual CouplingSignal driveSignal() const override
    {
        return CouplingSignal::VM;
    }

    virtual void evaluateHotPathStateForCell
    (
        const label cellI,
        const scalar modelTime,
        const scalar driveSignal,
        const scalar lambda,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& algebraicValues
    ) const override;

    virtual bool rushLarsenParametersForCell
    (
        const label cellI,
        const label stateI,
        const scalarUList& stateValues,
        const scalarUList& rateValues,
        const scalarUList& algebraicValues,
        scalar& steadyState,
        scalar& tau
    ) const override;

    virtual scalar activeTensionFromHotPath
    (
        const label cellI,
        const scalarUList& stateValues,
        const scalarUList& algebraicValues
    ) const override
    {
        return stateValues[::Ta];
    }

public:

    OverrideTypeName("GoktepeKuhlBatched");

    GoktepeKuhlBatched(const dictionary& dict, const label num);

    virtual ~GoktepeKuhlBatched() = default;

    virtual Requirements requirements() const override
    {
        Requirements r;
        r.needVm = true;
        return r;
    }
};

} // End namespace Foam

#endif
