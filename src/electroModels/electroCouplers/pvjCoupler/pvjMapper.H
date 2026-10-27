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
    pvjMapper

Description
    Mapper for PVJ coupling. Each terminal exchanges with the myocardium
    cells whose centre lies within radius R of it, weighted by a radial
    kernel w(r) and normalised by the discrete sum of w|K| over those cells:

        uniform   w = 1
        linear    w = 1 - r/R
        gaussian  w = exp(-4.5 r^2/R^2), truncated at R
        cosine    w = (1 + cos(pi r/R))/2

    The cosine kernel is the only one whose value and slope vanish at R.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef pvjMapper_H
#define pvjMapper_H

#include "volFields.H"

namespace Foam
{

class pvjMapper
{
    const fvMesh& mesh_;
    pointField terminalLocations_;
    scalar radius_;
    labelList terminalCellIDs_;
    List<labelList> terminalCellSets_;
    List<scalarList> terminalCellWeights_;
    scalarField sphereVolumes_;

public:

    //- Kernel names accepted by the constructor
    static const wordList kernelTypes;

    //- Kernel weight at distance r for support radius R (r <= R)
    static scalar kernelWeight(const word& kernelType, scalar r, scalar R);

    pvjMapper
    (
        const fvMesh& mesh,
        const pointField& terminalLocations,
        scalar radius,
        const word& kernelType = "uniform"
    );

    //- Local cells in the support of each terminal
    const List<labelList>& terminalCellSets() const
    {
        return terminalCellSets_;
    }

    //- Kernel weight of each local support cell
    const List<scalarList>& terminalCellWeights() const
    {
        return terminalCellWeights_;
    }

    //- Global sum of w|K| over each support: the discrete normalisation
    const scalarField& supportNormalisation() const
    {
        return sphereVolumes_;
    }

    void gatherVm3DPvjs
    (
        const volScalarField& Vm,
        scalarField& values
    ) const;

    //- Per terminal: earliest cell activation newer than newerThan (-1 if none) and latest cell activation.
    void gatherActivationTimes
    (
        const volScalarField& activationTimeField,
        const scalarField& newerThan,
        scalarField& terminalActivationTime,
        scalarField& latestActivationTime
    ) const;

    void depositCoupling
    (
        const scalarField& couplingCurrent,
        volScalarField& sourceField
    ) const;

    void depositImplicitCoupling
    (
        const scalarField& networkVm,
        const scalarField& resistance,
        volScalarField& sourceField,
        volScalarField& implicitSourceCoeff
    ) const;

    void depositActivationTimes
    (
        const scalarField& terminalActivationTime,
        volScalarField& activationTimeField
    ) const;

    void volumetricSource
    (
        const scalarField& couplingCurrent,
        scalarField& source
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
