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

Description
    Utilities for transmural heterogeneity weights
    (endo / M-cell / epi) for spatially varying ionic properties.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef ionicHeterogeneity_H
#define ionicHeterogeneity_H

#include "scalar.H"
#include "word.H"

namespace Foam
{
namespace ionicHeterogeneity
{

struct TransmuralBandWeights
{
    scalar endo;
    scalar mCell;
    scalar epi;
};

void validateTransmuralBandConfig
(
    scalar endoMInterface,
    scalar mEpiInterface,
    scalar transitionWidth,
    const word& smoothing
);

scalar smoothingWeight(scalar x, const word& smoothing);

TransmuralBandWeights transmuralBandWeights
(
    scalar t,
    scalar endoMInterface,
    scalar mEpiInterface,
    scalar transitionWidth,
    const word& smoothing
);

} // End namespace ionicHeterogeneity
} // End namespace Foam

#endif

// ************************************************************************* //
