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
    ecgElectrodeCells

Description
    Cell of each electrode point, owned by exactly one processor.

    The cell is the containing cell from findCell. A point outside the mesh,
    or found by more than one processor, takes the nearest cell centre; ties
    go to the lowest processor. Any other outcome is fatal.

SourceFiles
    ecgElectrodeCells.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef ecgElectrodeCells_H
#define ecgElectrodeCells_H

#include "fvMesh.H"
#include "labelList.H"

namespace Foam
{

//- Owner cell of each point on mesh, unique across processors (-1 where not
//  owned)
labelList ecgElectrodeCells(const fvMesh& mesh, const List<vector>& points);

} // End namespace Foam

#endif

// ************************************************************************* //
