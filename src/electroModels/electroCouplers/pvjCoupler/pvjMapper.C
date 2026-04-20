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

#include <cmath>

namespace Foam
{

PVJMapper::PVJMapper
(
    const fvMesh& mesh,
    const pointField& terminalLocations,
    scalar radius,
    bool reportSetup
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

        label nearestCellI = -1;
        scalar nearestDistance = GREAT;

        forAll(centres, cellI)
        {
            const scalar distance =
                mag(centres[cellI] - terminalLocations_[i]);

            if (distance < nearestDistance)
            {
                nearestDistance = distance;
                nearestCellI = cellI;
            }
        }

        scalar globalNearestDistance = nearestDistance;
        reduce(globalNearestDistance, minOp<scalar>());

        scalar nearestCellVolume = GREAT;
        if
        (
            nearestCellI >= 0
         && mag(nearestDistance - globalNearestDistance) <= SMALL
        )
        {
            nearestCellVolume = cellVolumes[nearestCellI];
        }
        reduce(nearestCellVolume, minOp<scalar>());

        if (nearestCellVolume >= GREAT/2.0)
        {
            FatalErrorInFunction
                << "Could not find a myocardium cell near PVJ terminal " << i
                << " at " << terminalLocations_[i] << "."
                << exit(FatalError);
        }

        const scalar nearestCellLength = std::cbrt(nearestCellVolume);
        if (radius_ < nearestCellLength)
        {
            FatalErrorInFunction
                << "pvjRadius is smaller than the local myocardium cell size "
                << "for PVJ terminal " << i << "." << nl
                << "  terminalLocation        = " << terminalLocations_[i] << nl
                << "  pvjRadius               = " << radius_ << nl
                << "  nearestCellDistance     = " << globalNearestDistance << nl
                << "  nearestCellVolume       = " << nearestCellVolume << nl
                << "  equivalentCellLength    = " << nearestCellLength << nl
                << "Increase pvjRadius to at least the local cell length, "
                << "or refine the myocardium mesh near the PVJ."
                << exit(FatalError);
        }

        if (terminalCellIDs_[i] < 0)
        {
            terminalCellIDs_[i] =
            (
                nearestCellI >= 0
             && mag(nearestDistance - globalNearestDistance) <= SMALL
              ? nearestCellI
              : -1
            );
        }

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
        label nLocalCells = terminalCellSets_[i].size();
        label nGlobalCells = nLocalCells;
        reduce(nGlobalCells, sumOp<label>());

        reduce(sphereVolumes_[i], sumOp<scalar>());

        if (sphereVolumes_[i] <= SMALL)
        {
            sphereVolumes_[i] = 0.0;

            if
            (
                nearestCellI >= 0
             && mag(nearestDistance - globalNearestDistance) <= SMALL
            )
            {
                terminalCellIDs_[i] = nearestCellI;
                terminalCellSets_[i].setSize(1);
                terminalCellSets_[i][0] = nearestCellI;
                sphereVolumes_[i] = cellVolumes[nearestCellI];
            }
            else
            {
                terminalCellIDs_[i] = -1;
                terminalCellSets_[i].clear();
            }

            reduce(sphereVolumes_[i], sumOp<scalar>());
            nLocalCells = terminalCellSets_[i].size();
            nGlobalCells = nLocalCells;
            reduce(nGlobalCells, sumOp<label>());

            if (sphereVolumes_[i] <= SMALL)
            {
                FatalErrorInFunction
                    << "PVJ sphere volume is zero for junction " << i << ". "
                    << "Increase pvjRadius or move the PVJ near the mesh."
                    << exit(FatalError);
            }
        }

        if (reportSetup && Pstream::master())
        {
            Info<< "PVJ mapping " << i
                << ": location=" << terminalLocations_[i]
                << ", pvjRadius=" << radius_
                << ", nearestCellDistance=" << globalNearestDistance
                << ", equivalentCellLength=" << nearestCellLength
                << ", mappedCells=" << nGlobalCells
                << ", mappedVolume=" << sphereVolumes_[i]
                << ", sourcePerAmp=" << 1.0/sphereVolumes_[i]
                << nl;
        }
    }

    if (reportSetup && Pstream::master())
    {
        Info<< endl;
    }
}


void PVJMapper::gatherVm3DPvjs
(
    const volScalarField& Vm,
    scalarField& values
) const
{
    values.setSize(terminalLocations_.size());
    values = -GREAT;

    forAll(values, i)
    {
        if (terminalCellIDs_[i] >= 0)
        {
            values[i] = Vm[terminalCellIDs_[i]];
        }

        reduce(values[i], maxOp<scalar>());

        if (values[i] <= -GREAT/2.0)
        {
            values[i] = 0.0;
        }
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


void PVJMapper::depositActivationTimes
(
    const scalarField& terminalActivationTime,
    volScalarField& activationTimeField
) const
{
    if (terminalActivationTime.size() != terminalCellSets_.size())
    {
        FatalErrorInFunction
            << "Expected " << terminalCellSets_.size()
            << " PVJ activation-time values but received "
            << terminalActivationTime.size()
            << exit(FatalError);
    }

    scalarField& activationValues = activationTimeField.primitiveFieldRef();

    forAll(terminalCellSets_, i)
    {
        const scalar terminalTime = terminalActivationTime[i];

        if (terminalTime < 0.0)
        {
            continue;
        }

        forAll(terminalCellSets_[i], localI)
        {
            const label cellI = terminalCellSets_[i][localI];

            if
            (
                activationValues[cellI] < 0.0
             || terminalTime < activationValues[cellI]
            )
            {
                activationValues[cellI] = terminalTime;
            }
        }
    }

    activationTimeField.correctBoundaryConditions();
}


void PVJMapper::gatherActivationTimes
(
    const volScalarField& activationTimeField,
    scalarField& terminalActivationTime
) const
{
    terminalActivationTime.setSize(terminalCellSets_.size());
    terminalActivationTime = GREAT;

    const scalarField& activationValues = activationTimeField.primitiveField();

    forAll(terminalCellSets_, i)
    {
        forAll(terminalCellSets_[i], localI)
        {
            const label cellI = terminalCellSets_[i][localI];
            const scalar t = activationValues[cellI];

            if (t >= 0.0 && t < terminalActivationTime[i])
            {
                terminalActivationTime[i] = t;
            }
        }

        reduce(terminalActivationTime[i], minOp<scalar>());

        if (terminalActivationTime[i] >= GREAT/2.0)
        {
            terminalActivationTime[i] = -1.0;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
