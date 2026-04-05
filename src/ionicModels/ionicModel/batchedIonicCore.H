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
    Foam::batchedIonicCore

Description
    Compact storage core for batched GPU-oriented ionic models. It owns the
    SoA state/rate/algebraic/support buffers and the cell gather/scatter
    helpers used by the host-side execution wrapper.
\*---------------------------------------------------------------------------*/

#ifndef batchedIonicCore_H
#define batchedIonicCore_H

#include "scalarField.H"

namespace Foam
{

class batchedIonicCore
{
    label nCells_;
    label nStates_;
    label nAlgebraics_;
    label nHotPathSupport_;

    mutable scalarField statesSoA_;
    mutable scalarField ratesSoA_;
    mutable scalarField algebraicsSoA_;
    mutable scalarField supportSoA_;

    inline label offset(const label cellI, const label varI) const
    {
        return varI*nCells_ + cellI;
    }

public:

    batchedIonicCore
    (
        const label nCells,
        const label nStates,
        const label nAlgebraics
    )
    :
        nCells_(nCells),
        nStates_(nStates),
        nAlgebraics_(nAlgebraics),
        nHotPathSupport_(0),
        statesSoA_(nStates*nCells, 0.0),
        ratesSoA_(nStates*nCells, 0.0),
        algebraicsSoA_(nAlgebraics*nCells, 0.0),
        supportSoA_()
    {}

    inline label nCells() const
    {
        return nCells_;
    }

    inline label nStates() const
    {
        return nStates_;
    }

    inline label nAlgebraics() const
    {
        return nAlgebraics_;
    }

    inline label nHotPathSupport() const
    {
        return nHotPathSupport_;
    }

    void setAlgebraicsStorageEnabled(const bool enabled)
    {
        if (enabled)
        {
            algebraicsSoA_.setSize(nAlgebraics_*nCells_);
            algebraicsSoA_ = 0.0;
        }
        else
        {
            algebraicsSoA_.clear();
        }
    }

    void setHotPathSupportSize(const label nSupport)
    {
        nHotPathSupport_ = max(nSupport, label(0));
        supportSoA_.setSize(nHotPathSupport_*nCells_);
        supportSoA_ = 0.0;
    }

    inline scalar& stateRef(const label cellI, const label stateI) const
    {
        return statesSoA_[offset(cellI, stateI)];
    }

    inline scalar& rateRef(const label cellI, const label stateI) const
    {
        return ratesSoA_[offset(cellI, stateI)];
    }

    inline scalar& algebraicRef(const label cellI, const label algI) const
    {
        return algebraicsSoA_[offset(cellI, algI)];
    }

    inline scalar& supportRef(const label cellI, const label supportI) const
    {
        return supportSoA_[offset(cellI, supportI)];
    }

    void gatherCellState
    (
        const label cellI,
        scalarField& stateValues
    ) const
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            stateValues[stateI] = stateRef(cellI, stateI);
        }
    }

    void scatterCellState
    (
        const label cellI,
        const scalarField& stateValues
    ) const
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            stateRef(cellI, stateI) = stateValues[stateI];
        }
    }

    void scatterCellRatesAndAlgebraics
    (
        const label cellI,
        const scalarUList& rateValues,
        const scalarUList& algebraicValues,
        const bool persistAlgebraics
    ) const
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            rateRef(cellI, stateI) = rateValues[stateI];
        }

        if (!persistAlgebraics)
        {
            return;
        }

        for (label algI = 0; algI < nAlgebraics_; ++algI)
        {
            algebraicRef(cellI, algI) = algebraicValues[algI];
        }
    }

    void scatterCellRates
    (
        const label cellI,
        const scalarUList& rateValues
    ) const
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            rateRef(cellI, stateI) = rateValues[stateI];
        }
    }

    void scatterCellSupport
    (
        const label cellI,
        const scalarUList& supportValues
    ) const
    {
        for (label supportI = 0; supportI < supportValues.size(); ++supportI)
        {
            supportRef(cellI, supportI) = supportValues[supportI];
        }
    }

    void clearTransientSolveData(const bool persistAlgebraics) const
    {
        ratesSoA_ = 0.0;
        supportSoA_ = 0.0;

        if (persistAlgebraics)
        {
            algebraicsSoA_ = 0.0;
        }
    }

    scalarField& statesField() const
    {
        return statesSoA_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
