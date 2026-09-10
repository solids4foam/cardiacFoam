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
    Mapper for PVJ coupling.

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

    pvjMapper
    (
        const fvMesh& mesh,
        const pointField& terminalLocations,
        scalar radius,
        const word& kernelType = "uniform"
    );

    void gatherVm3DPvjs
    (
        const volScalarField& Vm,
        scalarField& values
    ) const;

    void gatherActivationTimes
    (
        const volScalarField& activationTimeField,
        scalarField& terminalActivationTime
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
