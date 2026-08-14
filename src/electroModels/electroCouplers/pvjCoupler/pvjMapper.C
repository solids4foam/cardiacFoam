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

\*---------------------------------------------------------------------------*/

#include "pvjMapper.H"

#include "PstreamReduceOps.H"

#include <cmath>

namespace Foam
{

pvjMapper::pvjMapper
(
    const fvMesh& mesh,
    const pointField& terminalLocations,
    scalar radius,
    const word& kernelType
)
:
    mesh_(mesh),
    terminalLocations_(terminalLocations),
    radius_(radius),
    terminalCellIDs_(terminalLocations.size(), -1),
    terminalCellSets_(terminalLocations.size()),
    terminalCellWeights_(terminalLocations.size()),
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
        DynamicList<scalar> cellsInRadiusWeights;
        forAll(centres, cellI)
        {
            const scalar r = mag(centres[cellI] - terminalLocations_[i]);
            if (r <= radius_)
            {
                scalar w = 1.0;
                if (kernelType == "gaussian")
                {
                    w = std::exp(-4.5 * r * r / (radius_ * radius_ + VSMALL));
                }
                else if (kernelType == "linear")
                {
                    w = 1.0 - r / (radius_ + VSMALL);
                }

                cellsInRadius.append(cellI);
                cellsInRadiusWeights.append(w);
                sphereVolumes_[i] += w * cellVolumes[cellI];
            }
        }

        terminalCellSets_[i] = cellsInRadius;
        terminalCellWeights_[i] = cellsInRadiusWeights;
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
                terminalCellWeights_[i].setSize(1);
                terminalCellWeights_[i][0] = 1.0;
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

    }
}


void pvjMapper::gatherVm3DPvjs
(
    const volScalarField& Vm,
    scalarField& values
) const
{
    values.setSize(terminalLocations_.size());
    values = 0.0;

    forAll(values, i)
    {
        scalar localWeightedSum = 0.0;
        forAll(terminalCellSets_[i], localI)
        {
            const label cellI = terminalCellSets_[i][localI];
            localWeightedSum += Vm[cellI] * mesh_.V()[cellI] * terminalCellWeights_[i][localI];
        }

        values[i] = localWeightedSum;
        reduce(values[i], sumOp<scalar>());

        if (sphereVolumes_[i] > SMALL)
        {
            values[i] /= sphereVolumes_[i];
        }
        else
        {
            values[i] = 0.0;
        }
    }
}


void pvjMapper::volumetricSource
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


void pvjMapper::depositCoupling
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
            source[terminalCellSets_[i][localI]] += volumetric[i] * terminalCellWeights_[i][localI];
        }
    }

    sourceField.correctBoundaryConditions();
}


void pvjMapper::depositImplicitCoupling
(
    const scalarField& networkVm,
    const scalarField& resistance,
    volScalarField& sourceField,
    volScalarField& implicitSourceCoeff
) const
{
    if
    (
        networkVm.size() != sphereVolumes_.size()
     || resistance.size() != sphereVolumes_.size()
    )
    {
        FatalErrorInFunction
            << "Expected " << sphereVolumes_.size()
            << " PVJ implicit coupling values but received "
            << networkVm.size() << " voltages and "
            << resistance.size() << " resistances."
            << exit(FatalError);
    }

    scalarField& source = sourceField.primitiveFieldRef();
    scalarField& coeff = implicitSourceCoeff.primitiveFieldRef();

    forAll(terminalCellSets_, i)
    {
        const scalar sourcePerVoltage =
            1.0/(resistance[i]*sphereVolumes_[i]);

        forAll(terminalCellSets_[i], localI)
        {
            const label cellI = terminalCellSets_[i][localI];
            const scalar weight = terminalCellWeights_[i][localI];

            source[cellI] += networkVm[i]*sourcePerVoltage*weight;
            coeff[cellI] += sourcePerVoltage*weight;
        }
    }

    sourceField.correctBoundaryConditions();
    implicitSourceCoeff.correctBoundaryConditions();
}


void pvjMapper::depositActivationTimes
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


void pvjMapper::gatherActivationTimes
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
