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

#include "ionicHeterogeneity.H"
#include "error.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ionicHeterogeneity::validateTransmuralBandConfig
(
    const scalar endoMInterface,
    const scalar mEpiInterface,
    const scalar transitionWidth,
    const word& smoothing,
    const word& transitionMode
)
{
    if (transitionMode != "blend" && transitionMode != "hard")
    {
        FatalErrorInFunction
            << "Unsupported ionicHeterogeneity transitionMode '"
            << transitionMode
            << "'. Supported transitionMode values: blend, hard."
            << exit(FatalError);
    }

    if (smoothing != "smoothstep")
    {
        FatalErrorInFunction
            << "Unsupported ionicHeterogeneity smoothing '" << smoothing
            << "'. Supported smoothing: smoothstep."
            << exit(FatalError);
    }

    if
    (
        endoMInterface <= 0.0
     || mEpiInterface >= 1.0
     || endoMInterface >= mEpiInterface
     || transitionWidth < 0.0
     || (
            transitionMode == "blend"
         && (
                endoMInterface + transitionWidth > mEpiInterface + SMALL
             || mEpiInterface + transitionWidth > 1.0 + SMALL
            )
        )
    )
    {
        FatalErrorInFunction
            << "Invalid transmuralBands configuration: "
            << "endoMInterface=" << endoMInterface
            << ", mEpiInterface=" << mEpiInterface
            << ", transitionWidth=" << transitionWidth
            << ", transitionMode=" << transitionMode
            << ". Expected 0 < endoMInterface < mEpiInterface < 1 "
            << "and, for transitionMode blend, non-overlapping "
            << "transition bands."
            << exit(FatalError);
    }
}


Foam::scalar Foam::ionicHeterogeneity::smoothingWeight
(
    const scalar x,
    const word& smoothing
)
{
    if (smoothing == "smoothstep")
    {
        return x*x*(3.0 - 2.0*x);
    }

    FatalErrorInFunction
        << "Unsupported ionicHeterogeneity smoothing '" << smoothing
        << "'. Supported smoothing: smoothstep."
        << exit(FatalError);

    return 0.0;
}


Foam::ionicHeterogeneity::TransmuralBandWeights
Foam::ionicHeterogeneity::transmuralBandWeights
(
    const scalar t,
    const scalar endoMInterface,
    const scalar mEpiInterface,
    const scalar transitionWidth,
    const word& smoothing,
    const word& transitionMode
)
{
    if (transitionMode == "hard" || transitionWidth <= SMALL)
    {
        if (t <= endoMInterface)
        {
            return {1.0, 0.0, 0.0};
        }

        if (t <= mEpiInterface)
        {
            return {0.0, 1.0, 0.0};
        }

        return {0.0, 0.0, 1.0};
    }

    const scalar endoMUpper = endoMInterface + transitionWidth;
    const scalar mEpiUpper = mEpiInterface + transitionWidth;

    if (t <= endoMInterface)
    {
        return {1.0, 0.0, 0.0};
    }

    if (t < endoMUpper)
    {
        const scalar x = (t - endoMInterface)/transitionWidth;
        const scalar w = smoothingWeight(x, smoothing);
        return {1.0 - w, w, 0.0};
    }

    if (t <= mEpiInterface)
    {
        return {0.0, 1.0, 0.0};
    }

    if (t < mEpiUpper)
    {
        const scalar x = (t - mEpiInterface)/transitionWidth;
        const scalar w = smoothingWeight(x, smoothing);
        return {0.0, 1.0 - w, w};
    }

    return {0.0, 0.0, 1.0};
}

// ************************************************************************* //
