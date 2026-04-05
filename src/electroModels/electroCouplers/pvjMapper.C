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

#include "pvjMapper.H"

#include "PstreamReduceOps.H"

namespace Foam
{

PVJMapper::PVJMapper
(
    const fvMesh& mesh,
    const pointField& terminalLocations,
    scalar radius
)
:
    mesh_(mesh),
    terminalLocations_(terminalLocations),
    radius_(radius),
    terminalCellIDs_(terminalLocations.size(), -1),
    terminalCellSets_(terminalLocations.size()),
    sphereVolumes_(terminalLocations.size(), 0.0)
{
    const vectorField& centres = mesh_.C();
    const scalarField& cellVolumes = mesh_.V();

    forAll(terminalLocations_, i)
    {
        terminalCellIDs_[i] = mesh_.findCell(terminalLocations_[i]);

        DynamicList<label> cellsInRadius;
        forAll(centres, cellI)
        {
            if (mag(centres[cellI] - terminalLocations_[i]) <= radius_)
            {
                cellsInRadius.append(cellI);
                sphereVolumes_[i] += cellVolumes[cellI];
            }
        }

        terminalCellSets_[i] = cellsInRadius;
        reduce(sphereVolumes_[i], sumOp<scalar>());

        if (sphereVolumes_[i] <= SMALL)
        {
            FatalErrorInFunction
                << "PVJ sphere volume is zero for junction " << i << ". "
                << "Increase pvjRadius or move the PVJ inside the mesh."
                << exit(FatalError);
        }
    }
}


void PVJMapper::gatherVm3DPvjs
(
    const volScalarField& Vm,
    scalarField& values
) const
{
    values.setSize(terminalLocations_.size());
    values = 0.0;

    forAll(values, i)
    {
        if (terminalCellIDs_[i] >= 0)
        {
            values[i] = Vm[terminalCellIDs_[i]];
        }

        reduce(values[i], maxOp<scalar>());
    }
}


void PVJMapper::volumetricSource
(
    const scalarField& couplingCurrent,
    scalarField& source
) const
{
    if (couplingCurrent.size() != sphereVolumes_.size())
    {
        FatalErrorInFunction
            << "Expected " << sphereVolumes_.size()
            << " PVJ coupling values but received " << couplingCurrent.size()
            << exit(FatalError);
    }

    source.setSize(couplingCurrent.size());
    forAll(source, i)
    {
        source[i] = couplingCurrent[i]/sphereVolumes_[i];
    }
}


void PVJMapper::depositCoupling
(
    const scalarField& couplingCurrent,
    volScalarField& sourceField
) const
{
    scalarField volumetric;
    volumetricSource(couplingCurrent, volumetric);
    scalarField& source = sourceField.primitiveFieldRef();

    forAll(terminalCellSets_, i)
    {
        forAll(terminalCellSets_[i], localI)
        {
            source[terminalCellSets_[i][localI]] += volumetric[i];
        }
    }

    sourceField.correctBoundaryConditions();
}

} // End namespace Foam

// ************************************************************************* //
