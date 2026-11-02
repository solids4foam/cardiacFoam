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

labelList ecgElectrodeCells
(
    const fvMesh& mesh,
    const wordList& names,
    const List<vector>& points
)
{
    labelList cells(points.size(), -1);

    forAll(points, pointI)
    {
        label cellI = mesh.findCell(points[pointI]);
        const label ownerCount =
            returnReduce(label(cellI >= 0 ? 1 : 0), sumOp<label>());

        if (ownerCount == 0)
        {
            FatalErrorInFunction
                << "Electrode '" << names[pointI] << "' at " << points[pointI]
                << " lies outside mesh '" << mesh.name() << "': no cell "
                << "contains it. Torso ECG electrodes must lie in the mesh."
                << exit(FatalError);
        }

        if (ownerCount > 1)
        {
            const scalar distanceSqr =
                cellI >= 0
              ? Foam::magSqr(mesh.C().primitiveField()[cellI] - points[pointI])
              : GREAT;
            const scalar nearestSqr =
                returnReduce(distanceSqr, minOp<scalar>());

            const bool isNearest
            (
                cellI >= 0
             && Foam::mag(distanceSqr - nearestSqr)
              <= max(SMALL, ROOTSMALL*max(scalar(1), nearestSqr))
            );
            const label ownerProc = returnReduce
            (
                isNearest ? Pstream::myProcNo() : Pstream::nProcs(),
                minOp<label>()
            );

            if (Pstream::myProcNo() != ownerProc)
            {
                cellI = -1;
            }

            Info<< "ecgElectrodeCells: electrode '" << names[pointI]
                << "' at " << points[pointI] << " lies in a cell of "
                << ownerCount << " processors; using the cell of processor "
                << ownerProc << " at distance " << Foam::sqrt(nearestSqr)
                << " on mesh '" << mesh.name() << "'." << nl;
        }

        cells[pointI] = cellI;
    }

    return cells;
}

} // End namespace Foam

// ************************************************************************* //
