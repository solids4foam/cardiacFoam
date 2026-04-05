/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef PVJMapper_H
#define PVJMapper_H

#include "volFields.H"

namespace Foam
{

class PVJMapper
{
    const fvMesh& mesh_;
    pointField terminalLocations_;
    scalar radius_;
    labelList terminalCellIDs_;
    List<labelList> terminalCellSets_;
    scalarField sphereVolumes_;

public:

    PVJMapper
    (
        const fvMesh& mesh,
        const pointField& terminalLocations,
        scalar radius
    );

    void gatherVm3DPvjs
    (
        const volScalarField& Vm,
        scalarField& values
    ) const;

    void depositCoupling
    (
        const scalarField& couplingCurrent,
        volScalarField& sourceField
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
