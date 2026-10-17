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

#include "ecgElectrodeCells.H"
#include "volFields.H"
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


labelList ecgElectrodeCells(const fvMesh& mesh, const List<vector>& points)
{
    labelList cells(points.size(), -1);

    forAll(points, pointI)
    {
        cells[pointI] = mesh.findCell(points[pointI]);

        label ownerCount = cells[pointI] >= 0 ? 1 : 0;
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
                        mesh,
                        points[pointI],
                        localDistanceSqr
                    );
            }
            else if (cells[pointI] >= 0)
            {
                localCandidate = cells[pointI];
                localDistanceSqr =
                    Foam::magSqr
                    (
                        mesh.C().primitiveField()[localCandidate]
                      - points[pointI]
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

            cells[pointI] =
            (
                isNearest
             && Pstream::myProcNo() == ownerProc
              ? localCandidate
              : -1
            );

            ownerCount = cells[pointI] >= 0 ? 1 : 0;
            reduce(ownerCount, sumOp<label>());

            if (cells[pointI] >= 0)
            {
                if (initialOwnerCount == 0)
                {
                    Info<< "ecgElectrodeCells: electrode at " << points[pointI]
                        << " is outside findCell ownership; using nearest cell "
                        << cells[pointI]
                        << " at distance " << Foam::sqrt(globalDistanceSqr)
                        << " on mesh '" << mesh.name() << "'."
                        << nl;
                }
                else
                {
                    Info<< "ecgElectrodeCells: electrode at " << points[pointI]
                        << " has ambiguous parallel findCell ownership; "
                        << "using nearest owner cell "
                        << cells[pointI]
                        << " at distance " << Foam::sqrt(globalDistanceSqr)
                        << " on mesh '" << mesh.name() << "'."
                        << nl;
                }
            }
        }

        if (ownerCount != 1)
        {
            FatalErrorInFunction
                << "Electrode at " << points[pointI]
                << " must be owned by exactly one processor/cell in mesh '"
                << mesh.name() << "', but owner count is "
                << ownerCount << "."
                << exit(FatalError);
        }
    }

    return cells;
}

} // End namespace Foam

// ************************************************************************* //
