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
    gpuMath

Description
    Small branch-free helpers that are friendlier to batched CPU/GPU ionic
    kernels. These are intentionally local utilities and do not replace the
    exact OpenFOAM-style `heaviside` helper used elsewhere in the code base.
\*---------------------------------------------------------------------------*/

#ifndef gpuMath_H
#define gpuMath_H

#include "scalar.H"

#include <cmath>

namespace Foam
{

//- Smooth approximation of the Heaviside step function.
//  The default width keeps the transition narrow while avoiding a hard branch.
inline scalar smoothHeaviside(const scalar x, const scalar width = scalar(1.0e-4))
{
    return scalar(0.5) * (scalar(1.0) + std::tanh(x / width));
}

//- Smooth approximation of a thresholded indicator function.
inline scalar smoothHeaviside
(
    const scalar x,
    const scalar threshold,
    const scalar width
)
{
    return smoothHeaviside(x - threshold, width);
}

} // End namespace Foam

#endif

// ************************************************************************* //
