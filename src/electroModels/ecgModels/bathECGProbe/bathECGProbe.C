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

#include "bathECGProbe.H"

#include "ecgDomain.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

namespace Foam
{

namespace
{

label nearestCell
(
    const fvMesh& mesh,
    const point& samplePoint,
    scalar& nearestDistanceSqr
)
{
    label nearestCellI = -1;
    nearestDistanceSqr = GREAT;

    const vectorField& centres = mesh.C().primitiveField();

    forAll(centres, cellI)
    {
        const scalar distSqr = Foam::magSqr(centres[cellI] - samplePoint);

        if (distSqr < nearestDistanceSqr)
        {
            nearestDistanceSqr = distSqr;
            nearestCellI = cellI;
        }
    }

    return nearestCellI;
}

} // End anonymous namespace

defineTypeNameAndDebug(bathECGProbe, 0);
addToRunTimeSelectionTable(ecgSolver, bathECGProbe, dictionary);


bathECGProbe::bathECGProbe(const dictionary& dict)
:
    reportElectrodeLookup_
    (
        dict.lookupOrDefault<Switch>("reportElectrodeLookup", true)
    ),
    electrodeCells_(),
    cellsBuilt_(false)
{}


void bathECGProbe::buildElectrodeCells(const ecgDomain& domain)
{
    const volScalarField* phiEPtr = domain.phiEPtr();

    if (!phiEPtr)
    {
        FatalErrorInFunction
            << "bathECGProbe requires a global phiE field from the ECG "
            << "state provider."
            << exit(FatalError);
    }

    const fvMesh& phiEMesh = phiEPtr->mesh();
    const List<vector>& electrodes = domain.electrodePositions();

    electrodeCells_.setSize(electrodes.size(), -1);

    forAll(electrodes, electrodeI)
    {
        electrodeCells_[electrodeI] = phiEMesh.findCell(electrodes[electrodeI]);

        label ownerCount = electrodeCells_[electrodeI] >= 0 ? 1 : 0;
        reduce(ownerCount, sumOp<label>());

        if (ownerCount == 0)
        {
            scalar localDistanceSqr = GREAT;
            const label localNearest =
                nearestCell(phiEMesh, electrodes[electrodeI], localDistanceSqr);

            scalar globalDistanceSqr = localDistanceSqr;
            reduce(globalDistanceSqr, minOp<scalar>());

            electrodeCells_[electrodeI] =
            (
                localNearest >= 0
             && Foam::mag(localDistanceSqr - globalDistanceSqr)
              <= max(SMALL, ROOTSMALL*max(scalar(1), globalDistanceSqr))
              ? localNearest
              : -1
            );

            ownerCount = electrodeCells_[electrodeI] >= 0 ? 1 : 0;
            reduce(ownerCount, sumOp<label>());

            if (reportElectrodeLookup_ && electrodeCells_[electrodeI] >= 0)
            {
                Info<< "bathECGProbe: electrode "
                    << domain.electrodeNames()[electrodeI]
                    << " at " << electrodes[electrodeI]
                    << " is outside findCell ownership; using nearest cell "
                    << electrodeCells_[electrodeI]
                    << " at distance " << Foam::sqrt(globalDistanceSqr)
                    << " on phiE mesh '" << phiEMesh.name() << "'."
                    << nl;
            }
        }

        if (ownerCount != 1)
        {
            FatalErrorInFunction
                << "Electrode " << domain.electrodeNames()[electrodeI]
                << " at " << electrodes[electrodeI]
                << " must be owned by exactly one processor/cell in phiE "
                << "mesh '" << phiEMesh.name() << "', but owner count is "
                << ownerCount << "."
                << exit(FatalError);
        }
    }

    if (reportElectrodeLookup_)
    {
        Info<< "bathECGProbe: located " << electrodes.size()
            << " electrodes on phiE mesh '" << phiEMesh.name() << "'."
            << endl;
    }

    cellsBuilt_ = true;
}


void bathECGProbe::solve
(
    ecgDomain& domain,
    scalar t0,
    scalar dt,
    scalarField& values
)
{
    (void)t0;
    (void)dt;

    if (!cellsBuilt_)
    {
        buildElectrodeCells(domain);
    }

    const volScalarField* phiEPtr = domain.phiEPtr();

    if (!phiEPtr)
    {
        FatalErrorInFunction
            << "bathECGProbe requires phiE, but domain.phiEPtr() returned null."
            << exit(FatalError);
    }

    const scalarField& phiE = phiEPtr->primitiveField();

    values.setSize(electrodeCells_.size());
    values = 0.0;

    forAll(electrodeCells_, electrodeI)
    {
        const label cellI = electrodeCells_[electrodeI];
        if (cellI >= 0)
        {
            values[electrodeI] = phiE[cellI];
        }
    }

    if (Pstream::parRun())
    {
        Pstream::listCombineGather(values, plusEqOp<scalar>());
        Pstream::listCombineScatter(values);
    }
}

} // End namespace Foam

// ************************************************************************* //
