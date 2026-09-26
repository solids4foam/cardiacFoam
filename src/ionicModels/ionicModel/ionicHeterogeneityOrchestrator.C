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
    Spatial ionic heterogeneity configuration functions for ionic models.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "ionicHeterogeneityOrchestrator.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "ionicSelector.H"
#include "error.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ionicHeterogeneityOrchestrator::configureGradientAxisHeterogeneity
(
    const ionicModel& model,
    const scalarField& fieldValues,
    const dictionary& dict,
    PtrList<scalarField>& heterogeneousConstants
)
{
    const word axisName(dict.dictName());

    const wordList requiredKeys({"beta", "scalingMin", "scalingMax", "variables"});
    for (const word& key : requiredKeys)
    {
        if (!dict.found(key))
        {
            FatalErrorInFunction
                << "gradientAxes '" << axisName << "' for ionic model "
                << model.type() << " has no '" << key << "' entry. "
                << "Required: field, beta, scalingMin, scalingMax, variables."
                << exit(FatalError);
        }
    }

    const scalar beta(dict.get<scalar>("beta"));
    const scalar scalingMin(dict.get<scalar>("scalingMin"));
    const scalar scalingMax(dict.get<scalar>("scalingMax"));
    const wordList variables(dict.lookup("variables"));

    ionicHeterogeneity::validateGradientAxisConfig
    (
        axisName, scalingMin, scalingMax, variables, model.type()
    );

    const label nConst = model.ioNumConstants();
    const char* const* names = model.ioConstantNames();
    const scalarField* baseConstants = model.ioConstantsPtr();

    if (!names || nConst <= 0 || !baseConstants || baseConstants->empty())
    {
        FatalErrorInFunction
            << "gradientAxes '" << axisName << "' was requested for ionic "
            << "model " << model.type() << ", but this model does not "
            << "expose constant metadata (ioConstantNames / ioConstantsPtr)."
            << exit(FatalError);
    }

    labelList indices(variables.size(), -1);
    forAll(variables, vi)
    {
        for (label ci = 0; ci < nConst; ci++)
        {
            if (word(names[ci]) == variables[vi])
            {
                indices[vi] = ci;
                break;
            }
        }
        if (indices[vi] < 0)
        {
            FatalErrorInFunction
                << "gradientAxes '" << axisName << "': variable '"
                << variables[vi] << "' not found in constant names of "
                << "ionic model " << model.type() << ". Available constants: ";
            for (label ci = 0; ci < nConst; ci++)
            {
                FatalErrorInFunction << names[ci] << ' ';
            }
            FatalErrorInFunction << exit(FatalError);
        }
    }

    if (heterogeneousConstants.empty())
    {
        heterogeneousConstants.setSize(fieldValues.size());
        forAll(fieldValues, cellI)
        {
            heterogeneousConstants.set(cellI, new scalarField(*baseConstants));
        }
    }
    else if (heterogeneousConstants.size() != fieldValues.size())
    {
        FatalErrorInFunction
            << "gradientAxes '" << axisName << "': field has "
            << fieldValues.size() << " values but " << model.type()
            << " has " << heterogeneousConstants.size()
            << " heterogeneous constant sets."
            << exit(FatalError);
    }

    forAll(fieldValues, cellI)
    {
        const scalar d = min(max(fieldValues[cellI], scalar(0.0)), scalar(1.0));
        const scalar f =
            ionicHeterogeneity::apexBaseScale(d, beta, scalingMin, scalingMax);

        scalarField& consts = heterogeneousConstants[cellI];
        forAll(indices, vi)
        {
            consts[indices[vi]] *= f;
        }
    }
}


void Foam::ionicHeterogeneityOrchestrator::configureRegionHeterogeneity
(
    const ionicModel& model,
    const scalarField& transmuralDistance,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
)
{
    const auto* statesPtr = model.ioStatesPtr();

    if (statesPtr && transmuralDistance.size() != statesPtr->size())
    {
        FatalErrorInFunction
            << "Transmural distance field has " << transmuralDistance.size()
            << " values, but " << model.type() << " was configured with "
            << statesPtr->size() << " integration points."
            << exit(FatalError);
    }

    if (!heterogeneityDict.found("mode"))
    {
        FatalErrorInFunction
            << "ionicHeterogeneity for ionic model " << model.type()
            << " has no 'mode' entry. 'mode' is required: namedRegions or "
            << "cellZoneRegions."
            << exit(FatalError);
    }

    const word mode(heterogeneityDict.lookup("mode"));

    if (mode == "namedRegions")
    {
        configureNamedRegionHeterogeneity
        (
            model, transmuralDistance, heterogeneityDict, heterogeneousConstants,
            heterogeneousInitialStates
        );
        return;
    }

    if (mode == "cellZoneRegions")
    {
        configureCellZoneRegionHeterogeneity
        (
            model, transmuralDistance, heterogeneityDict, heterogeneousConstants,
            heterogeneousInitialStates
        );
        return;
    }

    FatalErrorInFunction
        << "Unsupported " << model.type() << " ionicHeterogeneity mode '"
        << mode << "'. Supported modes: namedRegions, cellZoneRegions."
        << exit(FatalError);
}


void Foam::ionicHeterogeneityOrchestrator::configureNamedRegionHeterogeneity
(
    const ionicModel& model,
    const scalarField& fieldValues,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
)
{
    if (!heterogeneityDict.found("transitionMode"))
    {
        FatalErrorInFunction
            << "ionicHeterogeneity mode namedRegions requires a "
            << "'transitionMode' entry for ionic model " << model.type()
            << ". Supported: blend, hard."
            << exit(FatalError);
    }

    const word transitionMode(heterogeneityDict.lookup("transitionMode"));

    if (transitionMode != "blend" && transitionMode != "hard")
    {
        FatalErrorInFunction
            << "Unsupported ionicHeterogeneity transitionMode '"
            << transitionMode << "' for mode namedRegions. Supported: "
            << "blend, hard."
            << exit(FatalError);
    }

    word smoothing;
    scalar transitionWidth = 0.0;

    if (transitionMode == "blend")
    {
        if (!heterogeneityDict.found("transitionWidth"))
        {
            FatalErrorInFunction
                << "ionicHeterogeneity mode namedRegions with "
                << "transitionMode blend requires a 'transitionWidth' "
                << "entry for ionic model " << model.type() << "."
                << exit(FatalError);
        }

        if (!heterogeneityDict.found("smoothing"))
        {
            FatalErrorInFunction
                << "ionicHeterogeneity mode namedRegions with "
                << "transitionMode blend requires a 'smoothing' entry "
                << "for ionic model " << model.type() << "."
                << exit(FatalError);
        }

        transitionWidth = heterogeneityDict.get<scalar>("transitionWidth");
        smoothing = word(heterogeneityDict.lookup("smoothing"));

        if (smoothing != "smoothstep")
        {
            FatalErrorInFunction
                << "Unsupported ionicHeterogeneity smoothing '" << smoothing
                << "' for mode namedRegions. Supported: smoothstep."
                << exit(FatalError);
        }

        if (transitionWidth < 0.0)
        {
            FatalErrorInFunction
                << "Invalid ionicHeterogeneity transitionWidth "
                << transitionWidth << " for mode namedRegions. Expected a "
                << "non-negative value."
                << exit(FatalError);
        }
    }

    if (!heterogeneityDict.found("regions"))
    {
        FatalErrorInFunction
            << "ionicHeterogeneity mode namedRegions requires a 'regions' "
            << "sub-dictionary for ionic model " << model.type() << "."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedFieldRegion> regions =
        ionicHeterogeneity::parseNamedFieldRegions
        (
            heterogeneityDict.subDict("regions")
        );

    if (transitionMode == "blend")
    {
        for (label regionI = 1; regionI < regions.size(); ++regionI)
        {
            const scalar regionWidth =
                regions[regionI].rangeMax - regions[regionI].rangeMin;

            if (transitionWidth > regionWidth + SMALL)
            {
                FatalErrorInFunction
                    << "Invalid ionicHeterogeneity transitionWidth "
                    << transitionWidth << " for mode namedRegions. Region "
                    << regions[regionI].name << " has width "
                    << regionWidth << ", so this transition would overlap "
                    << "the next boundary."
                    << exit(FatalError);
            }
        }
    }

    blendNamedRegions
    (
        model,
        fieldValues,
        regions,
        transitionWidth,
        smoothing,
        transitionMode,
        heterogeneousConstants,
        heterogeneousInitialStates
    );
}


void Foam::ionicHeterogeneityOrchestrator::blendNamedRegions
(
    const ionicModel& model,
    const scalarField& fieldValues,
    const List<ionicHeterogeneity::NamedFieldRegion>& regions,
    const scalar transitionWidth,
    const word& smoothing,
    const word& transitionMode,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
)
{
    wordList regionNames(regions.size());
    forAll(regions, i)
    {
        regionNames[i] = regions[i].name;
    }

    PtrList<scalarField> regionConstants(regions.size());
    PtrList<scalarField> regionStates(regions.size());
    bool blendStates = (heterogeneousInitialStates != nullptr);

    forAll(regions, i)
    {
        regionConstants.set
        (
            i,
            new scalarField
            (
                constantsForRegion(model, regionNames[i], regions[i].baseline, regionNames)
            )
        );

        if (regionConstants[i].empty())
        {
            FatalErrorInFunction
                << "ionicHeterogeneity region-based heterogeneity was "
                << "requested for ionic model " << model.type() << ", but "
                << "constantsForTissue() returned no constants. This model "
                << "does not support region-based heterogeneity."
                << exit(FatalError);
        }

        if (blendStates)
        {
            regionStates.set
            (
                i,
                new scalarField
                (
                    initialStatesForRegion
                    (
                        model, regionNames[i], regions[i].baseline, regionNames
                    )
                )
            );

            if (regionStates[i].empty())
            {
                blendStates = false;
            }
        }
    }

    heterogeneousConstants.clear();
    heterogeneousConstants.setSize(fieldValues.size());

    if (blendStates)
    {
        heterogeneousInitialStates->clear();
        heterogeneousInitialStates->setSize(fieldValues.size());
    }

    forAll(fieldValues, cellI)
    {
        const scalar rawT = fieldValues[cellI];

        if (rawT < -SMALL || rawT > 1.0 + SMALL)
        {
            FatalErrorInFunction
                << "Named-region field value t=" << rawT
                << " at integration point " << cellI
                << " is outside the expected [0, 1] range."
                << exit(FatalError);
        }

        const scalar t = min(max(rawT, scalar(0.0)), scalar(1.0));

        const List<ionicHeterogeneity::NamedRegionWeight> weights =
            ionicHeterogeneity::namedRegionWeightsAt
            (
                t, regions, transitionWidth, smoothing, transitionMode
            );

        scalarField mappedConstants(regionConstants[0].size(), 0.0);
        scalarField mappedStates;
        if (blendStates)
        {
            mappedStates.setSize(regionStates[0].size(), 0.0);
        }

        forAll(weights, wI)
        {
            const label regionI = regionNames.find(weights[wI].name);
            mappedConstants += weights[wI].weight*regionConstants[regionI];

            if (blendStates)
            {
                mappedStates += weights[wI].weight*regionStates[regionI];
            }
        }

        heterogeneousConstants.set(cellI, new scalarField(mappedConstants));

        if (blendStates)
        {
            heterogeneousInitialStates->set
            (
                cellI, new scalarField(mappedStates)
            );
        }
    }
}


Foam::scalarField Foam::ionicHeterogeneityOrchestrator::constantsForRegion
(
    const ionicModel& model,
    const word& regionName,
    const word& baseline,
    const wordList& knownRegionNames
)
{
    const wordList supported(model.supportedTissueTypes());
    if (!supported.found(baseline))
    {
        FatalErrorInFunction
            << "ionicHeterogeneity region '" << regionName
            << "' uses baseline '" << baseline << "', which ionic model "
            << model.type() << " does not support. Supported tissue types: "
            << supported
            << exit(FatalError);
    }

    const label tissueFlag = ionicSelector::tissueFlag(baseline);

    scalarField constants = model.constantsForTissue(tissueFlag);

    if (constants.empty())
    {
        return constants;
    }

    ionicModelIO::applyConstantOverrides
    (
        constants,
        model.ioConstantNames(),
        model.ioNumConstants(),
        model.dict(),
        model.type(),
        regionName,
        baseline,
        knownRegionNames
    );

    return constants;
}


Foam::scalarField Foam::ionicHeterogeneityOrchestrator::initialStatesForRegion
(
    const ionicModel& model,
    const word& regionName,
    const word& baseline,
    const wordList& knownRegionNames
)
{
    (void)regionName;
    (void)knownRegionNames;

    const label tissueFlag = ionicSelector::tissueFlag(baseline);

    return model.initialStatesForTissue(tissueFlag);
}


void Foam::ionicHeterogeneityOrchestrator::configureCellZoneRegionHeterogeneity
(
    const ionicModel& model,
    const scalarField& regionIndices,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
)
{
    const dictionary& regionsDict = heterogeneityDict.subDict("regions");
    const auto regions = ionicHeterogeneity::parseNamedCellZoneRegions(regionsDict);
    const label nRegions = regions.size();

    // Setup arrays of constant and initial state fields per region
    PtrList<scalarField> regionConstants(nRegions);
    PtrList<scalarField> regionInitialStates;
    if (heterogeneousInitialStates)
    {
        regionInitialStates.setSize(nRegions);
    }

    wordList regionNames(nRegions);
    forAll(regions, i)
    {
        regionNames[i] = regions[i].name;
    }

    forAll(regions, i)
    {
        const word& regionName = regions[i].name;
        const word& baseline = regions[i].baseline;

        regionConstants.set
        (
            i,
            new scalarField
            (
                constantsForRegion(model, regionName, baseline, regionNames)
            )
        );

        if (heterogeneousInitialStates)
        {
            regionInitialStates.set
            (
                i,
                new scalarField
                (
                    initialStatesForRegion(model, regionName, baseline, regionNames)
                )
            );
        }
    }

    const label nCells = regionIndices.size();

    heterogeneousConstants.clear();
    heterogeneousConstants.setSize(nCells);

    if (heterogeneousInitialStates)
    {
        heterogeneousInitialStates->clear();
        heterogeneousInitialStates->setSize(nCells);
    }

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        const label rIdx = round(regionIndices[cellI]);

        if (rIdx < 0 || rIdx >= nRegions)
        {
            FatalErrorInFunction
                << "Cell " << cellI << " mapped to region index " << rIdx
                << " but there are only " << nRegions << " cellZone regions defined."
                << exit(FatalError);
        }

        heterogeneousConstants.set(cellI, new scalarField(regionConstants[rIdx]));

        if (heterogeneousInitialStates)
        {
            heterogeneousInitialStates->set
            (
                cellI, new scalarField(regionInitialStates[rIdx])
            );
        }
    }
}

// ************************************************************************* //
