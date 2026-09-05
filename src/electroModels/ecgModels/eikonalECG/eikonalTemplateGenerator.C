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

#include "eikonalTemplateGenerator.H"
#include "ionicModel.H"
#include "ionicHeterogeneity.H"
#include "stimulusIO.H"

#include <cmath>

namespace Foam
{
namespace eikonalECG_templates
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


//- tPoints/anchorNames for mode transmuralBands: unchanged from the
//  original 3-anchor implementation.
void transmuralBandAnchors
(
    const dictionary& heterogeneityDict,
    scalarField& tPoints,
    wordList& anchorNames
)
{
    const scalar endoMInterface =
        heterogeneityDict.lookupOrDefault<scalar>("endoMInterface", 0.3);
    const scalar mEpiInterface =
        heterogeneityDict.lookupOrDefault<scalar>("mEpiInterface", 0.7);
    const scalar transitionWidth =
        heterogeneityDict.lookupOrDefault<scalar>("transitionWidth", 0.1);
    const word smoothing =
        heterogeneityDict.lookupOrDefault<word>("smoothing", "smoothstep");
    const word transitionMode =
        heterogeneityDict.lookupOrDefault<word>("transitionMode", "blend");

    ionicHeterogeneity::validateTransmuralBandConfig
    (
        endoMInterface,
        mEpiInterface,
        transitionWidth,
        smoothing,
        transitionMode
    );

    scalar tMid = 0.5*(endoMInterface + mEpiInterface);

    if (transitionMode == "blend" && transitionWidth > SMALL)
    {
        const scalar endoMUpper = endoMInterface + transitionWidth;

        if (endoMUpper >= mEpiInterface - SMALL)
        {
            FatalErrorInFunction
                << "eikonalTemplateGenerator: transitionWidth ("
                << transitionWidth << ") leaves no pure mid-myocardium "
                << "region between endoMInterface+transitionWidth ("
                << endoMUpper << ") and mEpiInterface (" << mEpiInterface
                << "). Reduce transitionWidth or widen the mCells band."
                << exit(FatalError);
        }

        tMid = 0.5*(endoMUpper + mEpiInterface);
    }

    tPoints.setSize(3);
    tPoints[0] = 0.0;
    tPoints[1] = tMid;
    tPoints[2] = 1.0;

    anchorNames.setSize(3);
    anchorNames[0] = "endocardium";
    anchorNames[1] = "mid-myocardium";
    anchorNames[2] = "epicardium";

    checkUnitInterval(tPoints[0], anchorNames[0]);
    checkUnitInterval(tPoints[1], anchorNames[1]);
    checkUnitInterval(tPoints[2], anchorNames[2]);
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
            << "eikonalTemplateGenerator: ionicHeterogeneity mode "
            << "namedRegions requires a 'regions' sub-dictionary."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedFieldRegion> regions =
        ionicHeterogeneity::parseNamedFieldRegions
        (
            heterogeneityDict.subDict("regions")
        );

    const word transitionMode =
        heterogeneityDict.lookupOrDefault<word>("transitionMode", "blend");
    const scalar transitionWidth =
        heterogeneityDict.lookupOrDefault<scalar>("transitionWidth", 0.1);

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
                << "eikonalTemplateGenerator: transitionWidth ("
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
            << "eikonalTemplateGenerator: ionicHeterogeneity mode "
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
    const scalar dt
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

    const word mode =
        heterogeneityDict.lookupOrDefault<word>("mode", "transmuralBands");

    scalarField tPoints;
    wordList anchorNames;

    if (mode == "transmuralBands")
    {
        transmuralBandAnchors(heterogeneityDict, tPoints, anchorNames);
    }
    else if (mode == "namedRegions")
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
            << "eikonalTemplateGenerator supports ionicHeterogeneity mode "
            << "'transmuralBands', 'namedRegions', or 'cellZoneRegions' "
            << "only; got '" << mode << "'."
            << exit(FatalError);
    }

    const label nPoints = tPoints.size();

    if (!ionicModelConfig.found("singleCellStimulus"))
    {
        FatalErrorInFunction
            << "eikonalTemplateGenerator requires a singleCellStimulus "
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

    // initialDeltaT is milliseconds; solveODE() below takes seconds.
    autoPtr<ionicModel> modelPtr =
        ionicModel::New(modelDict, nPoints, dt*1000.0, true);
    ionicModel& model = modelPtr();

    model.configureIonicHeterogeneity
    (
        tPoints,
        modelDict.subDict("ionicHeterogeneity")
    );

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

    for (label p = 0; p < nPoints; ++p)
    {
        result[p].times.setSize(nSamples);
        result[p].valuesMv.setSize(nSamples);
        result[p].times[0] = 0.0;
        // mV, not V -- unit conversion is the caller's job.
        result[p].valuesMv[0] = model.signal(p, CouplingSignal::VM);
    }

    for (label k = 1; k < nSamples; ++k)
    {
        model.solveODE(t, dt, Vm, Im);
        t += dt;

        for (label p = 0; p < nPoints; ++p)
        {
            result[p].times[k] = scalar(k)*dt;
            result[p].valuesMv[k] = model.signal(p, CouplingSignal::VM);
        }
    }

    for (label p = 0; p < nPoints; ++p)
    {
        validateTemplate(result[p], anchorNames[p]);
    }

    return result;
}

} // End namespace eikonalECG_templates
} // End namespace Foam

// ************************************************************************* //
