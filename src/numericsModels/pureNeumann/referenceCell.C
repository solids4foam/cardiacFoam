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

#include "referenceCell.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::locateReferenceCell(const fvMesh& mesh, const point& p)
{
    const label refCell = mesh.findCell(p);
    const label ownerCount =
        returnReduce(label(refCell >= 0 ? 1 : 0), sumOp<label>());

    if (ownerCount != 1)
    {
        FatalErrorInFunction
            << "phiERefPoint " << p << " must lie in exactly one cell of "
            << "one processor of mesh '" << mesh.name() << "', but "
            << ownerCount << " processors reported a containing cell."
            << exit(FatalError);
    }

    return refCell;
}


// ************************************************************************* //
