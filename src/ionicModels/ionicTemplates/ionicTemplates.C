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

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "ionicTemplates.H"
#include "ionicModel.H"
#include "ionicHeterogeneity.H"
#include "stimulusIO.H"
#include "restartStateIO.H"
#include "dictionaryKeyCheck.H"
#include "IFstream.H"
#include "volFields.H"

#include <cmath>

namespace Foam
{
namespace ionicTemplates
{

namespace
{

void checkUnitInterval(const scalar t, const word& anchor)
{
    if (t < 0.0 || t > 1.0)
    {
        FatalErrorInFunction
            << "Transmural reference point for anchor '" << anchor
            << "' is " << t << ", outside [0, 1]."
            << exit(FatalError);
    }
}


void validateTemplate(const DynamicTemplate& tpl, const word& anchor)
{
    if (tpl.times.empty() || tpl.valuesMv.empty())
    {
        FatalErrorInFunction
            << "Generated an empty action-potential template for anchor '"
            << anchor << "'."
            << exit(FatalError);
    }

    forAll(tpl.valuesMv, i)
    {
        if (!std::isfinite(tpl.valuesMv[i]) || !std::isfinite(tpl.times[i]))
        {
            FatalErrorInFunction
                << "Non-finite sample in the generated action-potential "
                << "template for anchor '" << anchor << "' at index " << i
                << " (t=" << tpl.times[i] << ", Vm=" << tpl.valuesMv[i]
                << ")."
                << exit(FatalError);
        }

        if (i > 0 && tpl.times[i] <= tpl.times[i - 1])
        {
            FatalErrorInFunction
                << "Generated action-potential template for anchor '"
                << anchor << "' has non-increasing times at index " << i
                << " (" << tpl.times[i - 1] << " -> " << tpl.times[i] << ")."
                << exit(FatalError);
        }
    }

    const scalar resting = tpl.valuesMv[0];
    scalar peak = resting;
    forAll(tpl.valuesMv, i)
    {
        peak = max(peak, tpl.valuesMv[i]);
    }
    const scalar amplitude = peak - resting;

    if (amplitude < 50.0)
    {
        FatalErrorInFunction
            << "Generated action-potential template for anchor '" << anchor
            << "' has amplitude " << amplitude << " mV, below the required "
            << "50 mV minimum (resting=" << resting << " mV, peak=" << peak
            << " mV)."
            << exit(FatalError);
    }

    const scalar lastValue = tpl.valuesMv.last();
    const scalar recovery = mag(lastValue - resting);

    if (recovery > 5.0)
    {
        FatalErrorInFunction
            << "Generated action-potential template for anchor '" << anchor
            << "' did not recover to within 5 mV of its resting value by "
            << "capture end (resting=" << resting << " mV, final="
            << lastValue << " mV, |final-resting|=" << recovery << " mV)."
            << exit(FatalError);
    }
}


//- tPoints/anchorNames for mode namedRegions: one anchor per region, at
//  the midpoint of that region's "pure" sub-range (the part of its range
//  not eaten into by a blend zone with the PREVIOUS region -- mirroring
//  transmuralBandAnchors' own mid-band treatment, generalized to N
//  regions). Region 0 has no left neighbour, so its pure range is its
//  full range.
void namedRegionAnchors
(
    const dictionary& heterogeneityDict,
    scalarField& tPoints,
    wordList& anchorNames
)
{
    if (!heterogeneityDict.found("regions"))
    {
        FatalErrorInFunction
            << "ionicTemplates: ionicHeterogeneity mode "
            << "namedRegions requires a 'regions' sub-dictionary."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedFieldRegion> regions =
        ionicHeterogeneity::parseNamedFieldRegions
        (
            heterogeneityDict.subDict("regions")
        );

    if (!heterogeneityDict.found("transitionMode"))
    {
        FatalErrorInFunction
            << "ionicTemplates: ionicHeterogeneity mode "
            << "namedRegions requires a 'transitionMode' entry."
            << exit(FatalError);
    }

    const word transitionMode(heterogeneityDict.get<word>("transitionMode"));

    scalar transitionWidth = 0.0;

    if (transitionMode == "blend")
    {
        if (!heterogeneityDict.found("transitionWidth"))
        {
            FatalErrorInFunction
                << "ionicTemplates: ionicHeterogeneity mode "
                << "namedRegions with transitionMode blend requires a "
                << "'transitionWidth' entry."
                << exit(FatalError);
        }

        transitionWidth = heterogeneityDict.get<scalar>("transitionWidth");
    }
    else if (transitionMode != "hard")
    {
        FatalErrorInFunction
            << "Unsupported ionicHeterogeneity transitionMode '"
            << transitionMode << "'. Supported: blend, hard."
            << exit(FatalError);
    }

    const label nRegions = regions.size();
    tPoints.setSize(nRegions);
    anchorNames.setSize(nRegions);

    forAll(regions, i)
    {
        anchorNames[i] = regions[i].name;

        if (i == 0 || transitionMode != "blend" || transitionWidth <= SMALL)
        {
            tPoints[i] = 0.5*(regions[i].rangeMin + regions[i].rangeMax);
            continue;
        }

        const scalar pureStart = regions[i].rangeMin + transitionWidth;

        if (pureStart >= regions[i].rangeMax - SMALL)
        {
            FatalErrorInFunction
                << "ionicTemplates: transitionWidth ("
                << transitionWidth << ") leaves no pure region for '"
                << regions[i].name << "' (range " << regions[i].rangeMin
                << " " << regions[i].rangeMax << "). Reduce transitionWidth "
                << "or widen this region."
                << exit(FatalError);
        }

        tPoints[i] = 0.5*(pureStart + regions[i].rangeMax);
        checkUnitInterval(tPoints[i], anchorNames[i]);
    }
}

//- tPoints/anchorNames for mode cellZoneRegions: one anchor per region, no
//  blending (this mode has no transitionWidth/smoothing/transitionMode at
//  all in the monodomain path). Anchor i is given the synthetic coordinate
//  scalar(i); ionicHeterogeneityOrchestrator::configureCellZoneRegionHeterogeneity
//  resolves a cell's region via round(regionIndices[cellI]) used as a
//  direct 0-based array index, so scalar(i) correctly selects regions[i].
void cellZoneRegionAnchors
(
    const dictionary& heterogeneityDict,
    scalarField& tPoints,
    wordList& anchorNames
)
{
    if (!heterogeneityDict.found("regions"))
    {
        FatalErrorInFunction
            << "ionicTemplates: ionicHeterogeneity mode "
            << "cellZoneRegions requires a 'regions' sub-dictionary."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedCellZoneRegion> regions =
        ionicHeterogeneity::parseNamedCellZoneRegions
        (
            heterogeneityDict.subDict("regions")
        );

    const label nRegions = regions.size();
    tPoints.setSize(nRegions);
    anchorNames.setSize(nRegions);

    forAll(regions, i)
    {
        tPoints[i] = scalar(i);
        anchorNames[i] = regions[i].name;
    }
}

} // End unnamed namespace


List<DynamicTemplate> generatePersonalizedTemplates
(
    const dictionary& ionicModelConfig,
    const dictionary& heterogeneityDict,
    const label nBeats,
    const scalar captureDuration,
    const scalar dt,
    const dictionary* regionStates,
    const fileName& stateRoot
)
{
    if (nBeats < 1 || captureDuration <= 0.0 || dt <= 0.0)
    {
        FatalErrorInFunction
            << "Invalid template generation controls: nBeats=" << nBeats
            << ", captureDuration=" << captureDuration << ", dt=" << dt
            << ". Require nBeats >= 1 and captureDuration, dt > 0."
            << exit(FatalError);
    }

    if (!heterogeneityDict.found("mode"))
    {
        FatalErrorInFunction
            << "ionicTemplates requires ionicHeterogeneity 'mode' "
            << "to be set: 'namedRegions' or 'cellZoneRegions'."
            << exit(FatalError);
    }

    const word mode(heterogeneityDict.get<word>("mode"));

    scalarField tPoints;
    wordList anchorNames;

    if (mode == "namedRegions")
    {
        namedRegionAnchors(heterogeneityDict, tPoints, anchorNames);
    }
    else if (mode == "cellZoneRegions")
    {
        cellZoneRegionAnchors(heterogeneityDict, tPoints, anchorNames);
    }
    else
    {
        FatalErrorInFunction
            << "ionicTemplates supports ionicHeterogeneity mode "
            << "'namedRegions' or 'cellZoneRegions' only; got '" << mode
            << "'."
            << exit(FatalError);
    }

    const label nPoints = tPoints.size();

    if (!ionicModelConfig.found("singleCellStimulus"))
    {
        FatalErrorInFunction
            << "ionicTemplates requires a singleCellStimulus "
            << "sub-dictionary defining an S1 pacing protocol (stim_start, "
            << "stim_period_S1, stim_duration, stim_amplitude); none was "
            << "found in ionicModelConfig."
            << exit(FatalError);
    }

    dictionary modelDict(ionicModelConfig);
    modelDict.set("ionicHeterogeneity", heterogeneityDict);

    dictionary stimDict(modelDict.subDict("singleCellStimulus"));
    stimDict.set("nstim1", nBeats);
    modelDict.set("singleCellStimulus", stimDict);

    const StimulusProtocol stim = stimulusIO::loadStimulusProtocol(modelDict);
    const scalar tCapture =
        (stim.stimStart + scalar(nBeats - 1)*stim.stimPeriodS1)*1e-3;

    autoPtr<ionicModel> modelPtr =
        ionicModel::New(modelDict, nPoints, dt, true);
    ionicModel& model = modelPtr();

    model.configureIonicHeterogeneity
    (
        tPoints,
        modelDict.subDict("ionicHeterogeneity")
    );

    if (regionStates)
    {
        const word modelName(modelDict.get<word>("ionicModel"));
        label stride = -1;
        List<scalarField> states(nPoints);
        forAll(anchorNames, p)
        {
            if (!regionStates->found(anchorNames[p]))
            {
                FatalIOErrorInFunction(*regionStates)
                    << "regionStates has no state file for region '"
                    << anchorNames[p] << "'. Regions: " << anchorNames
                    << exit(FatalIOError);
            }
            fileName statePath(regionStates->get<fileName>(anchorNames[p]));
            statePath.expand();
            if (!statePath.isAbsolute())
            {
                statePath = stateRoot/statePath;
            }
            states[p] =
                restartStateIO::readSingleCellState(statePath, modelName, stride);
        }
        if (!model.setStates(states))
        {
            FatalErrorInFunction
                << "Ionic model " << modelName << " did not accept the region "
                << "states from regionStates (" << stride << " states per cell)."
                << exit(FatalError);
        }
    }

    scalarField Vm(nPoints, 0.0);
    scalarField Im(nPoints, 0.0);

    scalar t = 0.0;
    const label nFullPreSteps = label(std::floor(tCapture/dt + 1e-9));
    for (label s = 0; s < nFullPreSteps; ++s)
    {
        model.solveODE(t, dt, Vm, Im);
        t += dt;
    }

    // Land exactly on tCapture, or the waveform shifts by up to one step.
    const scalar remainder = tCapture - t;
    if (remainder > SMALL)
    {
        model.solveODE(t, remainder, Vm, Im);
        t = tCapture;
    }

    const label nSamples = label(std::ceil(captureDuration/dt - 1e-9)) + 1;

    List<DynamicTemplate> result(nPoints);
    const bool hasCai = model.hasSignal(CouplingSignal::CAI);

    for (label p = 0; p < nPoints; ++p)
    {
        result[p].times.setSize(nSamples);
        result[p].valuesMv.setSize(nSamples);
        result[p].valuesCai.setSize(hasCai ? nSamples : 0);
        result[p].times[0] = 0.0;
        // mV, not V -- unit conversion is the caller's job.
        result[p].valuesMv[0] = model.signal(p, CouplingSignal::VM);
        if (hasCai)
        {
            result[p].valuesCai[0] = model.signal(p, CouplingSignal::CAI);
        }
    }

    for (label k = 1; k < nSamples; ++k)
    {
        model.solveODE(t, dt, Vm, Im);
        t += dt;

        for (label p = 0; p < nPoints; ++p)
        {
            result[p].times[k] = scalar(k)*dt;
            result[p].valuesMv[k] = model.signal(p, CouplingSignal::VM);
            if (hasCai)
            {
                result[p].valuesCai[k] = model.signal(p, CouplingSignal::CAI);
            }
        }
    }

    // Time origin at the steepest upstroke
    for (label p = 0; p < nPoints; ++p)
    {
        const scalarField& v = result[p].valuesMv;
        label kUp = 0;
        for (label k = 1; k < v.size() - 1; ++k)
        {
            if (v[k + 1] - v[k] > v[kUp + 1] - v[kUp])
            {
                kUp = k;
            }
        }
        const scalar tUp =
            0.5*(result[p].times[kUp] + result[p].times[kUp + 1]);
        result[p].times -= tUp;
    }

    for (label p = 0; p < nPoints; ++p)
    {
        validateTemplate(result[p], anchorNames[p]);
    }

    return result;
}


void checkControls(const dictionary& templatesDict)
{
    const wordList keys
    ({
        "dt", "duration", "nBeats", "ionicModelConfig", "regionStates"
    });
    wordRes allowed(keys.size());
    forAll(keys, i)
    {
        allowed[i] = wordRe(keys[i]);
    }
    dictionaryKeyCheck::checkKeys(templatesDict, allowed);

    if (!templatesDict.found("ionicModelConfig"))
    {
        FatalIOErrorInFunction(templatesDict)
            << "personalizedTemplates requires an ionicModelConfig "
            << "sub-dictionary." << exit(FatalIOError);
    }

    const dictionary& ionicModelConfig =
        templatesDict.subDict("ionicModelConfig");

    if (!ionicModelConfig.found("ionicModel"))
    {
        FatalIOErrorInFunction(templatesDict)
            << "personalizedTemplates.ionicModelConfig requires an "
            << "'ionicModel' entry." << exit(FatalIOError);
    }

    if (!ionicModelConfig.found("singleCellStimulus"))
    {
        FatalIOErrorInFunction(templatesDict)
            << "personalizedTemplates.ionicModelConfig requires a "
            << "'singleCellStimulus' sub-dictionary defining the S1 pacing "
            << "protocol." << exit(FatalIOError);
    }

    const label nBeats = templateBeats(templatesDict);
    const scalar duration = templatesDict.get<scalar>("duration");
    const scalar generationDt = templatesDict.get<scalar>("dt");

    if (nBeats < 1)
    {
        FatalIOErrorInFunction(templatesDict)
            << "personalizedTemplates.nBeats must be >= 1; got " << nBeats
            << "." << exit(FatalIOError);
    }

    if (duration <= 0.0 || generationDt <= 0.0)
    {
        FatalIOErrorInFunction(templatesDict)
            << "personalizedTemplates.duration and .dt must both be "
            << "positive; got duration=" << duration << ", dt="
            << generationDt << "." << exit(FatalIOError);
    }

    const StimulusProtocol stim =
        stimulusIO::loadStimulusProtocol(ionicModelConfig);

    if (duration > 1e-3*stim.stimPeriodS1)
    {
        FatalIOErrorInFunction(templatesDict)
            << "personalizedTemplates.duration (" << duration
            << " s) must not exceed one S1 period ("
            << 1e-3*stim.stimPeriodS1 << " s, from ionicModelConfig."
            << "singleCellStimulus.stim_period_S1=" << stim.stimPeriodS1
            << " ms): capture must stay within the final S1 beat."
            << exit(FatalIOError);
    }

    if (stim.nStim2 != 0)
    {
        FatalIOErrorInFunction(templatesDict)
            << "personalizedTemplates does not support S2 pacing "
            << "(ionicModelConfig.singleCellStimulus.nstim2 must be 0); "
            << "capture must contain exactly one unambiguous final S1 "
            << "response." << exit(FatalIOError);
    }
}


label templateBeats(const dictionary& templatesDict)
{
    return
        templatesDict.found("regionStates")
      ? templatesDict.lookupOrDefault<label>("nBeats", 1)
      : templatesDict.get<label>("nBeats");
}


List<DynamicTemplate> generate
(
    const dictionary& templatesDict,
    const dictionary& heterogeneityDict,
    const Time& runTime
)
{
    checkControls(templatesDict);

    // Region states in the setIonicRestartStateDict format
    autoPtr<dictionary> regionStatesPtr;
    if (templatesDict.found("regionStates"))
    {
        fileName mapPath(templatesDict.get<fileName>("regionStates"));
        mapPath.expand();
        if (!mapPath.isAbsolute())
        {
            mapPath = runTime.globalPath()/mapPath;
        }
        IFstream is(mapPath);
        if (!is.good())
        {
            FatalIOErrorInFunction(templatesDict)
                << "Cannot read regionStates file " << mapPath
                << exit(FatalIOError);
        }
        const dictionary mapDict(is);
        regionStatesPtr.reset(new dictionary(mapDict.subDict("regionStates")));
    }

    return generatePersonalizedTemplates
    (
        templatesDict.subDict("ionicModelConfig"),
        heterogeneityDict,
        templateBeats(templatesDict),
        templatesDict.get<scalar>("duration"),
        templatesDict.get<scalar>("dt"),
        regionStatesPtr.get(),
        runTime.globalPath()
    );
}


const dictionary* findHeterogeneityDict(const dictionary& electroProperties)
{
    const dictionary* hetDictPtr =
        electroProperties.findDict("ionicHeterogeneity");

    if (!hetDictPtr && electroProperties.found("myocardiumSolver"))
    {
        const word myocardiumSolver
        (
            electroProperties.get<word>("myocardiumSolver")
        );
        const dictionary* solverDictPtr =
            electroProperties.findDict(myocardiumSolver + "Coeffs");

        if (solverDictPtr)
        {
            hetDictPtr = solverDictPtr->findDict("ionicHeterogeneity");
        }
    }

    return hetDictPtr;
}


List<scalarField> regionWeights
(
    const fvMesh& mesh,
    const dictionary& heterogeneityDict
)
{
    const dictionary& hetDict = heterogeneityDict;
    const word mode(hetDict.get<word>("mode"));
    List<scalarField> weights;

    if (mode == "cellZoneRegions")
    {
        const List<ionicHeterogeneity::NamedCellZoneRegion> czRegions =
            ionicHeterogeneity::parseNamedCellZoneRegions
            (
                hetDict.subDict("regions")
            );

        weights.setSize(czRegions.size());
        forAll(weights, i)
        {
            weights[i].setSize(mesh.nCells(), 0.0);
        }

        boolList claimed(mesh.nCells(), false);

        forAll(czRegions, regionIndex)
        {
            const word& zoneName = czRegions[regionIndex].cellZone;
            const label zoneId = mesh.cellZones().findZoneID(zoneName);

            if (zoneId < 0)
            {
                FatalErrorInFunction
                    << "cellZone '" << zoneName << "' (region '"
                    << czRegions[regionIndex].name << "') not found in "
                    << "this mesh."
                    << exit(FatalError);
            }

            const labelList& zoneCells = mesh.cellZones()[zoneId];

            forAll(zoneCells, i)
            {
                const label cellI = zoneCells[i];

                if (claimed[cellI])
                {
                    FatalErrorInFunction
                        << "Cell " << cellI << " is claimed by more than "
                        << "one cellZoneRegions entry."
                        << exit(FatalError);
                }

                claimed[cellI] = true;
                weights[regionIndex][cellI] = 1.0;
            }
        }

        forAll(claimed, cellI)
        {
            if (!claimed[cellI])
            {
                FatalErrorInFunction
                    << "Cell " << cellI << " is not claimed by any "
                    << "cellZoneRegions entry. Every mesh cell must belong "
                    << "to exactly one named cell zone."
                    << exit(FatalError);
            }
        }

        return weights;
    }

    if (mode != "namedRegions")
    {
        FatalErrorInFunction
            << "Region weights support ionicHeterogeneity mode "
            << "'namedRegions' or 'cellZoneRegions' only; got '" << mode
            << "'."
            << exit(FatalError);
    }

    if (!hetDict.found("transitionMode"))
    {
        FatalErrorInFunction
            << "ionicHeterogeneity mode namedRegions requires a "
            << "'transitionMode' entry. Supported: blend, hard."
            << exit(FatalError);
    }

    const word transitionMode(hetDict.get<word>("transitionMode"));

    word smoothing;
    scalar transitionWidth = 0.0;

    if (transitionMode == "blend")
    {
        if (!hetDict.found("transitionWidth") || !hetDict.found("smoothing"))
        {
            FatalErrorInFunction
                << "ionicHeterogeneity mode namedRegions with "
                << "transitionMode blend requires both 'transitionWidth' "
                << "and 'smoothing' entries."
                << exit(FatalError);
        }

        transitionWidth = hetDict.get<scalar>("transitionWidth");
        smoothing = hetDict.get<word>("smoothing");
    }
    else if (transitionMode != "hard")
    {
        FatalErrorInFunction
            << "Unsupported ionicHeterogeneity transitionMode '"
            << transitionMode << "'. Supported: blend, hard."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedFieldRegion> regions =
        ionicHeterogeneity::parseNamedFieldRegions
        (
            hetDict.subDict("regions")
        );

    if (!hetDict.found("field"))
    {
        FatalErrorInFunction
            << "ionicHeterogeneity mode namedRegions requires a "
            << "'field' entry."
            << exit(FatalError);
    }

    const word fieldName(hetDict.get<word>("field"));
    autoPtr<volScalarField> tReadPtr;
    const volScalarField* tPtr = mesh.cfindObject<volScalarField>(fieldName);

    if (!tPtr)
    {
        // Transmural-distance fields are static anatomical inputs stored in
        // the latest time directory that contains the field.
        const word fieldInstance =
            mesh.time().findInstance(mesh.dbDir(), fieldName);

        tReadPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    fieldInstance,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        tPtr = &tReadPtr();
    }

    const scalarField& tField = tPtr->primitiveField();

    const label nRegions = regions.size();
    weights.setSize(nRegions);
    forAll(weights, i)
    {
        weights[i].setSize(mesh.nCells(), 0.0);
    }

    // Small, fixed-size lookup built once, not per-cell.
    List<word> regionNames(nRegions);
    forAll(regions, i)
    {
        regionNames[i] = regions[i].name;
    }

    forAll(tField, cellI)
    {
        const List<ionicHeterogeneity::NamedRegionWeight> sparse =
            ionicHeterogeneity::namedRegionWeightsAt
            (
                tField[cellI], regions, transitionWidth, smoothing, transitionMode
            );

        forAll(sparse, s)
        {
            forAll(regionNames, i)
            {
                if (regionNames[i] == sparse[s].name)
                {
                    weights[i][cellI] = sparse[s].weight;
                    break;
                }
            }
        }
    }

    return weights;
}

} // End namespace ionicTemplates
} // End namespace Foam

// ************************************************************************* //
