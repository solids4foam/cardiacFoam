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
#include <cmath>
#include <algorithm>
#include "DynamicList.H"
#include "scalarList.H"

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


Foam::scalar Foam::ionicHeterogeneity::apexBaseScale
(
    const scalar d,
    const scalar beta,
    const scalar scalingMin,
    const scalar scalingMax
)
{
    return scalingMin*(1.0 + (scalingMax/scalingMin - 1.0)*std::exp(-beta*d));
}


Foam::List<Foam::ionicHeterogeneity::NamedFieldRegion>
Foam::ionicHeterogeneity::synthesizeTransmuralBandRegions
(
    const scalar endoMInterface,
    const scalar mEpiInterface
)
{
    List<NamedFieldRegion> regions(3);

    regions[0].name = "endocardialCells";
    regions[0].rangeMin = 0.0;
    regions[0].rangeMax = endoMInterface;
    regions[0].baseline = "endocardialCells";

    regions[1].name = "mCells";
    regions[1].rangeMin = endoMInterface;
    regions[1].rangeMax = mEpiInterface;
    regions[1].baseline = "mCells";

    regions[2].name = "epicardialCells";
    regions[2].rangeMin = mEpiInterface;
    regions[2].rangeMax = 1.0;
    regions[2].baseline = "epicardialCells";

    return regions;
}


Foam::List<Foam::ionicHeterogeneity::NamedFieldRegion>
Foam::ionicHeterogeneity::parseNamedFieldRegions
(
    const dictionary& regionsDict
)
{
    DynamicList<NamedFieldRegion> regions;

    forAllConstIter(dictionary, regionsDict, iter)
    {
        const word regionName(iter().keyword());

        if (!iter().isDict())
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions entry '" << regionName
                << "' must be a dictionary."
                << exit(FatalError);
        }

        if (regionName == "global")
        {
            FatalErrorInFunction
                << "ionicHeterogeneity region name 'global' is reserved for "
                << "ionicConstantOverrides.global and cannot be used as a "
                << "named region."
                << exit(FatalError);
        }

        const dictionary& regionDict = regionsDict.subDict(regionName);

        if (!regionDict.found("range"))
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions." << regionName
                << " has no 'range' entry. Field-based regions must "
                << "declare 'range (min max);'."
                << exit(FatalError);
        }

        const scalarList range(regionDict.lookup("range"));

        if (range.size() != 2 || range[0] >= range[1])
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions." << regionName
                << ".range must be exactly two values (min max) with "
                << "min < max, got " << range
                << exit(FatalError);
        }

        static const wordList anatomicalNames
        {
            "epicardialCells", "mCells", "endocardialCells"
        };

        const word baseline = regionDict.lookupOrDefault<word>
        (
            "baseline",
            anatomicalNames.found(regionName) ? regionName : word("myocyte")
        );

        if (!anatomicalNames.found(baseline) && baseline != "myocyte")
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions." << regionName
                << ".baseline '" << baseline << "' is not a supported "
                << "tissue baseline. Supported: epicardialCells, mCells, "
                << "endocardialCells, myocyte."
                << exit(FatalError);
        }

        NamedFieldRegion region;
        region.name = regionName;
        region.rangeMin = range[0];
        region.rangeMax = range[1];
        region.baseline = baseline;
        regions.append(region);
    }

    if (regions.size() < 2)
    {
        FatalErrorInFunction
            << "ionicHeterogeneity mode namedRegions requires at least "
            << "two entries under 'regions', found " << regions.size()
            << exit(FatalError);
    }

    std::sort
    (
        regions.begin(), regions.end(),
        [](const NamedFieldRegion& a, const NamedFieldRegion& b)
        {
            return a.rangeMin < b.rangeMin;
        }
    );

    for (label i = 1; i < regions.size(); ++i)
    {
        if (mag(regions[i].rangeMin - regions[i - 1].rangeMax) > SMALL)
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions '" << regions[i - 1].name
                << "' (range " << regions[i - 1].rangeMin << " "
                << regions[i - 1].rangeMax << ") and '" << regions[i].name
                << "' (range " << regions[i].rangeMin << " "
                << regions[i].rangeMax << ") are not adjacent. Named "
                << "field regions must tile [0,1] with no gaps or overlaps."
                << exit(FatalError);
        }
    }

    if (mag(regions[0].rangeMin - 0.0) > SMALL)
    {
        FatalErrorInFunction
            << "ionicHeterogeneity.regions: first region '" << regions[0].name
            << "' starts at " << regions[0].rangeMin << ", not 0. Named "
            << "field regions must tile [0,1] with no gaps at the edges."
            << exit(FatalError);
    }

    if (mag(regions.last().rangeMax - 1.0) > SMALL)
    {
        FatalErrorInFunction
            << "ionicHeterogeneity.regions: last region '"
            << regions.last().name << "' ends at " << regions.last().rangeMax
            << ", not 1. Named field regions must tile [0,1] with no gaps "
            << "at the edges."
            << exit(FatalError);
    }

    return List<NamedFieldRegion>(regions);
}


Foam::List<Foam::ionicHeterogeneity::NamedCellZoneRegion>
Foam::ionicHeterogeneity::parseNamedCellZoneRegions
(
    const dictionary& regionsDict
)
{
    DynamicList<NamedCellZoneRegion> regions;
    DynamicList<word> seenZones;

    forAllConstIter(dictionary, regionsDict, iter)
    {
        const word regionName(iter().keyword());

        if (!iter().isDict())
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions entry '" << regionName
                << "' must be a dictionary."
                << exit(FatalError);
        }

        if (regionName == "global")
        {
            FatalErrorInFunction
                << "ionicHeterogeneity region name 'global' is reserved for "
                << "ionicConstantOverrides.global and cannot be used as a "
                << "named region."
                << exit(FatalError);
        }

        const dictionary& regionDict = regionsDict.subDict(regionName);

        if (!regionDict.found("cellZone"))
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions." << regionName
                << " has no 'cellZone' entry. Cell-zone based regions must "
                << "declare 'cellZone <word>;'."
                << exit(FatalError);
        }

        const word cellZoneName(regionDict.lookup("cellZone"));

        if (seenZones.found(cellZoneName))
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions." << regionName
                << ": cellZone '" << cellZoneName << "' is already claimed "
                << "by another region. Each cellZone may be assigned to "
                << "exactly one region."
                << exit(FatalError);
        }
        seenZones.append(cellZoneName);

        static const wordList anatomicalNames
        {
            "epicardialCells", "mCells", "endocardialCells"
        };

        const word baseline = regionDict.lookupOrDefault<word>
        (
            "baseline",
            anatomicalNames.found(regionName) ? regionName : word("myocyte")
        );

        if (!anatomicalNames.found(baseline) && baseline != "myocyte")
        {
            FatalErrorInFunction
                << "ionicHeterogeneity.regions." << regionName
                << ".baseline '" << baseline << "' is not a supported "
                << "tissue baseline. Supported: epicardialCells, mCells, "
                << "endocardialCells, myocyte."
                << exit(FatalError);
        }

        NamedCellZoneRegion region;
        region.name = regionName;
        region.cellZone = cellZoneName;
        region.baseline = baseline;
        regions.append(region);
    }

    if (regions.size() < 2)
    {
        FatalErrorInFunction
            << "ionicHeterogeneity mode cellZoneRegions requires at least "
            << "two entries under 'regions', found " << regions.size()
            << exit(FatalError);
    }

    return List<NamedCellZoneRegion>(regions);
}


Foam::List<Foam::ionicHeterogeneity::NamedRegionWeight>
Foam::ionicHeterogeneity::namedRegionWeightsAt
(
    const scalar t,
    const List<NamedFieldRegion>& regions,
    const scalar transitionWidth,
    const word& smoothing,
    const word& transitionMode
)
{
    const label nRegions = regions.size();

    label idx = nRegions - 1;
    for (label i = 0; i < nRegions; ++i)
    {
        // Exact boundaries belong to the lower field region, matching the
        // legacy transmuralBandWeights() convention for hard assignment and
        // zero-width transitions.
        if (t <= regions[i].rangeMax || i == nRegions - 1)
        {
            idx = i;
            break;
        }
    }

    if (transitionMode == "hard" || transitionWidth <= SMALL)
    {
        return List<NamedRegionWeight>(1, {regions[idx].name, scalar(1.0)});
    }

    if (idx > 0)
    {
        const scalar boundary = regions[idx - 1].rangeMax;
        if (t < boundary + transitionWidth)
        {
            const scalar x = (t - boundary)/transitionWidth;
            const scalar w = smoothingWeight(x, smoothing);
            List<NamedRegionWeight> out(2);
            out[0] = {regions[idx - 1].name, 1.0 - w};
            out[1] = {regions[idx].name, w};
            return out;
        }
    }

    return List<NamedRegionWeight>(1, {regions[idx].name, scalar(1.0)});
}

// ************************************************************************* //
