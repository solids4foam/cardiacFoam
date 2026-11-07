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

Function
    locateReferenceCell

Description
    The cell holding the reference point of a pure-Neumann potential: the
    cell of mesh that contains p, on the one processor whose cell contains
    it (-1 on the others). Fatal unless exactly one processor has such a
    cell. Used by the heart-only bidomain phiE and by the heart+torso
    conductor.

SourceFiles
    referenceCell.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef referenceCell_H
#define referenceCell_H

#include "fvMesh.H"

namespace Foam
{

label locateReferenceCell(const fvMesh& mesh, const point& p);

} // End namespace Foam

#endif

// ************************************************************************* //
