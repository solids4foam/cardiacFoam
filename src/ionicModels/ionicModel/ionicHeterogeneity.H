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
    (endo / M-cell / epi) and apex-to-base gradient scaling
    for spatially varying ionic properties.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef ionicHeterogeneity_H
#define ionicHeterogeneity_H

#include "scalar.H"
#include "word.H"
#include "List.H"
#include "dictionary.H"

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
    const word& smoothing,
    const word& transitionMode
);

//- Validate a gradientAxes axis's scalingMin/scalingMax/variables.
void validateGradientAxisConfig
(
    const word& axisName,
    scalar scalingMin,
    scalar scalingMax,
    const wordList& variables,
    const word& modelType
);

struct NamedFieldRegion
{
    word name;
    scalar rangeMin;
    scalar rangeMax;

    //- Explicit tissue baseline this region resolves to via
    //  constantsForTissue()/initialStatesForTissue(). If the region name
    //  differs from the baseline, ionicConstantOverrides.<name> is layered
    //  on top. One of:
    //  epicardialCells, mCells, endocardialCells, myocyte. Defaults to the
    //  region's own name if it is one of the three anatomical names,
    //  otherwise defaults to myocyte — set explicitly to override this.
    word baseline;
};

struct NamedCellZoneRegion
{
    word name;
    word cellZone;
    word baseline;
};

struct NamedRegionWeight
{
    word name;
    scalar weight;
};

//- Parse an open dictionary of {name: {range (min max); baseline <word>;}}
//  entries into a list sorted by rangeMin. 'baseline' is optional per
//  region (see NamedFieldRegion::baseline for the default rule) but if
//  given must be one of epicardialCells/mCells/endocardialCells/myocyte.
//  Fatals if any region is named 'global', if any region is missing
//  'range', if a range is malformed, if fewer than two regions are
//  declared, if 'baseline' is not a supported tissue name, or if the
//  sorted regions do not tile [0,1] contiguously with no gaps/overlaps
//  — including at the edges: the first region's rangeMin must be exactly
//  0 and the last region's rangeMax must be exactly 1 (within SMALL).
List<NamedFieldRegion> parseNamedFieldRegions(const dictionary& regionsDict);

//- Parse an open dictionary of {name: {cellZone <word>; baseline <word>;}}
//  entries. 'baseline' follows the same fallback rules as NamedFieldRegion.
//  Fatals if any region is missing 'cellZone'.
List<NamedCellZoneRegion> parseNamedCellZoneRegions(const dictionary& regionsDict);

//- Per-cell region weights (sum to 1) for field value t against a
//  pre-sorted, contiguous region list. transitionWidth <= SMALL or
//  transitionMode "hard" gives a single region with weight 1.
List<NamedRegionWeight> namedRegionWeightsAt
(
    scalar t,
    const List<NamedFieldRegion>& regions,
    scalar transitionWidth,
    const word& smoothing,
    const word& transitionMode
);

scalar smoothingWeight(scalar x, const word& smoothing);

TransmuralBandWeights transmuralBandWeights
(
    scalar t,
    scalar endoMInterface,
    scalar mEpiInterface,
    scalar transitionWidth,
    const word& smoothing,
    const word& transitionMode
);

scalar apexBaseScale
(
    scalar d,
    scalar beta,
    scalar scalingMin,
    scalar scalingMax
);

//- Expand a transmuralBands (endoMInterface, mEpiInterface) configuration
//  into the same three named regions
//  (endocardialCells, mCells, epicardialCells) that mode namedRegions
//  consumes directly, so callers that need a uniform List<NamedFieldRegion>
//  view of either mode (e.g. eikonalECG's personalizedTemplates) can treat
//  transmuralBands as the N=3 case of the same region-based machinery.
//  Mirrors the inline construction in
//  ionicHeterogeneityOrchestrator::configureTransmuralBandHeterogeneity();
//  callers must not assume this list is used to change that orchestrator's
//  own behavior, only to give a second consumer an equivalent view.
List<NamedFieldRegion> synthesizeTransmuralBandRegions
(
    scalar endoMInterface,
    scalar mEpiInterface
);

} // End namespace ionicHeterogeneity
} // End namespace Foam

#endif

// ************************************************************************* //
