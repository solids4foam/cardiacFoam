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
    Cell of each electrode of a torso ECG, owned by exactly one processor.

    The cell is the containing cell from findCell. An electrode outside the
    mesh is fatal and named. A point found by more than one processor (on a
    processor boundary) takes the containing cell whose centre is nearest;
    ties go to the lowest processor.

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

//- Owner cell of each electrode point on mesh, unique across processors
//  (-1 where not owned); names name the electrodes in the fatal error
labelList ecgElectrodeCells
(
    const fvMesh& mesh,
    const wordList& names,
    const List<vector>& points
);

} // End namespace Foam

#endif

// ************************************************************************* //
