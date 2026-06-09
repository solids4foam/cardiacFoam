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
    Foam::LandNiedererBatched

Description
    Land-Niederer 2017 active tension model implemented for batched GPU
    execution.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef LandNiedererBatched_H
#define LandNiedererBatched_H

#include "batchedActiveTensionModel.H"
#include "../LandNiederer/LandNiederer_2017Names.H"
#include "overrideTypeName.H"

namespace Foam
{

class LandNiedererBatched
:
    public batchedActiveTensionModel
{
    LandNiedererBatched(const LandNiedererBatched&) = delete;
    void operator=(const LandNiedererBatched&) = delete;

    mutable scalarField CONSTANTS_;
    mutable scalarField prevLambda_;
    mutable scalarField lambdaRate_;

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
        return CouplingSignal::CAI;
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

    virtual scalar activeTensionFromHotPath
    (
        const label cellI,
        const scalarUList& stateValues,
        const scalarUList& algebraicValues
    ) const override
    {
        return algebraicValues[AV_Ta];
    }

public:

    OverrideTypeName("LandNiedererBatched");

    LandNiedererBatched(const dictionary& dict, const label num);

    virtual ~LandNiedererBatched() = default;

    virtual Requirements requirements() const override
    {
        Requirements r;
        r.needCai = true;
        r.needVm = false;
        return r;
    }

    virtual scalar timeScaleFactor() const override
    {
        return 1000.0; // scales dt from s to ms for the batched executor
    }

    virtual void calculateTension
    (
        const scalar t,
        const scalar dt,
        const scalarField& lambda,
        scalarField& Ta
    ) override;
};

} // End namespace Foam

#endif
