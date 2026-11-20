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
    reportCompatibilityImbalance

Description
    Discrete imbalance of a pure-Neumann potential equation: the sum of its
    source and of its non-coupled boundary coefficients, which must vanish
    for the singular system to have a solution. The equation is not changed.

    The continuum problem is compatible. For a conservative operator with
    insulated walls the imbalance is round-off; a non-zero value comes from
    a non-conservative source, a wall that is not insulated, or unbalanced
    surface currents. The reference cell of the solve absorbs it, so it is
    reported, relative to the sum of the magnitudes of the same terms, every
    time the equation is assembled.

SourceFiles
    compatibilityImbalance.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef compatibilityImbalance_H
#define compatibilityImbalance_H

#include "fvMatrices.H"

namespace Foam
{

//- Report the discrete imbalance of the pure-Neumann equation eqn and
//  return it
scalar reportCompatibilityImbalance(const fvScalarMatrix& eqn);

}

#endif

// ************************************************************************* //
