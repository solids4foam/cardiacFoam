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

#include "torsoECG.H"

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

defineTypeNameAndDebug(torsoECG, 0);
addToRunTimeSelectionTable(ecgSolver, torsoECG, dictionary);


torsoECG::torsoECG(const dictionary& dict)
:
    electrodeCells_(),
    cellsBuilt_(false)
{}


void torsoECG::buildElectrodeCells(const ecgDomain& domain)
{
    const volScalarField* phiEPtr = domain.phiEPtr();

    if (!phiEPtr)
    {
        FatalErrorInFunction
            << "torsoECG requires a global phiE field from the ECG "
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
        const label initialOwnerCount = ownerCount;

        if (ownerCount != 1)
        {
            scalar localDistanceSqr = GREAT;
            label localCandidate = -1;

            if (ownerCount == 0)
            {
                localCandidate =
                    nearestCell
                    (
                        phiEMesh,
                        electrodes[electrodeI],
                        localDistanceSqr
                    );
            }
            else if (electrodeCells_[electrodeI] >= 0)
            {
                localCandidate = electrodeCells_[electrodeI];
                localDistanceSqr =
                    Foam::magSqr
                    (
                        phiEMesh.C().primitiveField()[localCandidate]
                      - electrodes[electrodeI]
                    );
            }

            scalar globalDistanceSqr = localDistanceSqr;
            reduce(globalDistanceSqr, minOp<scalar>());

            const bool isNearest
            (
                localCandidate >= 0
             && Foam::mag(localDistanceSqr - globalDistanceSqr)
              <= max(SMALL, ROOTSMALL*max(scalar(1), globalDistanceSqr))
            );

            label ownerProc = isNearest ? Pstream::myProcNo() : Pstream::nProcs();
            reduce(ownerProc, minOp<label>());

            electrodeCells_[electrodeI] =
            (
                isNearest
             && Pstream::myProcNo() == ownerProc
              ? localCandidate
              : -1
            );

            ownerCount = electrodeCells_[electrodeI] >= 0 ? 1 : 0;
            reduce(ownerCount, sumOp<label>());

            if (electrodeCells_[electrodeI] >= 0)
            {
                if (initialOwnerCount == 0)
                {
                    Info<< "torsoECG: electrode "
                        << domain.electrodeNames()[electrodeI]
                        << " at " << electrodes[electrodeI]
                        << " is outside findCell ownership; using nearest cell "
                        << electrodeCells_[electrodeI]
                        << " at distance " << Foam::sqrt(globalDistanceSqr)
                        << " on phiE mesh '" << phiEMesh.name() << "'."
                        << nl;
                }
                else
                {
                    Info<< "torsoECG: electrode "
                        << domain.electrodeNames()[electrodeI]
                        << " at " << electrodes[electrodeI]
                        << " has ambiguous parallel findCell ownership; "
                        << "using nearest owner cell "
                        << electrodeCells_[electrodeI]
                        << " at distance " << Foam::sqrt(globalDistanceSqr)
                        << " on phiE mesh '" << phiEMesh.name() << "'."
                        << nl;
                }
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

    Info<< "torsoECG: located " << electrodes.size()
        << " electrodes on phiE mesh '" << phiEMesh.name() << "'."
        << endl;

    cellsBuilt_ = true;
}


void torsoECG::solve
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
            << "torsoECG requires phiE, but domain.phiEPtr() returned null."
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
        Pstream::broadcast(values);
    }
}

} // End namespace Foam

// ************************************************************************* //
