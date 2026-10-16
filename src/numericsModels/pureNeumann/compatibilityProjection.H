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
    projectCompatibility

Description
    Removes the discrete imbalance of a pure-Neumann system before its
    reference cell is pinned: the sum of the source and of the non-coupled
    boundary coefficients is subtracted, weighted by cell volume.

SourceFiles
    compatibilityProjection.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef compatibilityProjection_H
#define compatibilityProjection_H

#include "fvMatrices.H"

namespace Foam
{

//- Remove the discrete imbalance of eqn and report it
void projectCompatibility(fvScalarMatrix& eqn);

}

#endif

// ************************************************************************* //
