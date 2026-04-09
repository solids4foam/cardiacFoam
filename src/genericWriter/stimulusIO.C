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

#include "stimulusIO.H"
#include "IOstreams.H"
#include "dimensionedScalar.H"
#include <cmath>

namespace Foam
{
namespace
{
const dictionary* singleCellStimulusDict(const dictionary& dict)
{
    if (!dict.found("singleCellStimulus"))
    {
        return nullptr;
    }

    return &dict.subDict("singleCellStimulus");
}

bool hasAnySingleCellStimulusKey(const dictionary& dict)
{
    return
        dict.found("stim_start")
     || dict.found("stim_period_S1")
     || dict.found("stim_duration")
     || dict.found("stim_amplitude")
     || dict.found("nstim1")
     || dict.found("stim_period_S2")
     || dict.found("nstim2");
}

const dictionary* externalStimulusDict(const dictionary& dict)
{
    if (!dict.found("externalStimulus"))
    {
        return nullptr;
    }

    return &dict.subDict("externalStimulus");
}

void readMonodomainStimulusBoxes
(
    const dictionary& stimDict,
    List<boundBox>& boxes
)
{
    const bool hasMinList = stimDict.found("stimulusLocationMinList");
    const bool hasMaxList = stimDict.found("stimulusLocationMaxList");

    if (hasMinList != hasMaxList)
    {
        FatalErrorInFunction
            << "Both stimulusLocationMinList and stimulusLocationMaxList "
            << "must be provided together."
            << abort(FatalError);
    }

    if (hasMinList)
    {
        const List<point> mins(stimDict.lookup("stimulusLocationMinList"));
        const List<point> maxs(stimDict.lookup("stimulusLocationMaxList"));

        if (mins.size() != maxs.size())
        {
            FatalErrorInFunction
                << "stimulusLocationMinList and stimulusLocationMaxList must "
                << "have the same size."
                << abort(FatalError);
        }

        boxes.setSize(mins.size());
        forAll(mins, i)
        {
            boxes[i] = boundBox(mins[i], maxs[i]);
        }
    }
    else
    {
        boxes.setSize(1);
        boxes[0] = boundBox
        (
            point(stimDict.lookup("stimulusLocationMin")),
            point(stimDict.lookup("stimulusLocationMax"))
        );
    }
}

void checkSize
(
    const List<scalar>& values,
    const label expected,
    const word& name
)
{
    if (values.size() != expected)
    {
        FatalErrorInFunction
            << name << " must have the same size as the stimulus box list."
            << abort(FatalError);
    }
}

void checkNonNegative
(
    const List<scalar>& values,
    const word& name
)
{
    forAll(values, i)
    {
        if (values[i] < 0.0)
        {
            FatalErrorInFunction
                << name << " must be non-negative. Found value "
                << values[i] << " at index " << i << "."
                << abort(FatalError);
        }
    }
}
}

StimulusProtocol stimulusIO::loadStimulusProtocol
(
    const dictionary& dict
)
{
    StimulusProtocol stim;

    const dictionary* stimDictPtr = singleCellStimulusDict(dict);
    if (!stimDictPtr)
    {
        // Default no-stimulus protocol if no explicit singleCellStimulus
        // dictionary is provided.
        return stim;
    }
    const dictionary& stimDict = *stimDictPtr;
    const bool hasAnyStimulusKey = hasAnySingleCellStimulusKey(stimDict);

    // Incomplete definitions are usually a typo and should fail fast.
    if
    (
        hasAnyStimulusKey
     && (
            !stimDict.found("stim_start")
         || !stimDict.found("stim_period_S1")
         || !stimDict.found("stim_duration")
         || !stimDict.found("stim_amplitude")
        )
    )
    {
        FatalErrorInFunction
            << "Incomplete single-cell stimulus protocol. Required keys are: "
            << "(stim_start, stim_period_S1, stim_duration, stim_amplitude)."
            << nl
            << "Define them inside sub-dictionary 'singleCellStimulus'."
            << exit(FatalError);
    }

    stim.stimStart = stimDict.lookupOrDefault<scalar>("stim_start", 0.0);
    stim.stimPeriodS1 = stimDict.lookupOrDefault<scalar>("stim_period_S1", 0.0);
    stim.stimDuration = stimDict.lookupOrDefault<scalar>("stim_duration", 0.0);
    stim.stimAmplitude = stimDict.lookupOrDefault<scalar>("stim_amplitude", 0.0);

    stim.nStim1 = stimDict.lookupOrDefault<label>
    (
        "nstim1",
        stim.stimPeriodS1 > SMALL ? 1 : 0
    );

    stim.stimPeriodS2 = stimDict.lookupOrDefault<scalar>("stim_period_S2", 0.0);
    stim.nStim2 = stimDict.lookupOrDefault<label>("nstim2", 0);

    return stim;
}

void stimulusIO::loadStimulusProtocol
(
    const dictionary& dict,
    scalarField& C,
    label stim_start,
    label stim_period_S1,
    label stim_duration,
    label stim_amplitude,
    label nstim1,
    label stim_period_S2,
    label nstim2
)
{
    const StimulusProtocol stim = loadStimulusProtocol(dict);
    C[stim_start] = stim.stimStart;
    C[stim_period_S1] = stim.stimPeriodS1;
    C[stim_duration] = stim.stimDuration;
    C[stim_amplitude] = stim.stimAmplitude;
    C[nstim1] = stim.nStim1;
    C[stim_period_S2] = stim.stimPeriodS2;
    C[nstim2] = stim.nStim2;
}

scalar stimulusIO::computeStimulus
(
    scalar VOI,
    const StimulusProtocol& stim
)
{
    return computeStimulus
    (
        VOI,
        stim.stimStart,
        stim.stimPeriodS1,
        stim.stimDuration,
        stim.stimAmplitude,
        stim.nStim1,
        stim.stimPeriodS2,
        stim.nStim2
    );
}

scalar stimulusIO::computeStimulus
(
    scalar VOI,
    scalar stim_start,
    scalar stim_period_S1,
    scalar stim_duration,
    scalar stim_amplitude,
    label  nstim1,
    scalar stim_period_S2,
    label  nstim2
)
{
    scalar Istim = 0.0;

    // ---------- S1 ----------
    if (stim_period_S1 > 0.0 && nstim1 > 0)
    {
        scalar tp = VOI - stim_start;

        if (tp >= 0 && tp <= stim_period_S1 * nstim1)
        {
            scalar phase = tp - floor(tp / stim_period_S1) * stim_period_S1;

            if (phase >= 0 && phase <= stim_duration)
            {
                Istim = -stim_amplitude;
            }
        }
    }

    // ---------- S2 ----------
    if (Istim == 0.0 && stim_period_S2 > 0.0 && nstim2 > 0)
    {
        scalar tS1End = stim_start + stim_period_S1 * nstim1;
        scalar tp = VOI - tS1End;

        if (tp >= 0 && tp <= stim_period_S2 * nstim2)
        {
            scalar phase = tp - floor(tp / stim_period_S2) * stim_period_S2;

            if (phase >= 0 && phase <= stim_duration)
            {
                Istim = -stim_amplitude;
            }
        }
    }

    return Istim;
}

bool stimulusIO::hasActiveStimulus(const StimulusProtocol& stim)
{
    const bool hasS1 =
        stim.stimPeriodS1 > SMALL
     && stim.nStim1 > 0
     && stim.stimDuration > SMALL
     && mag(stim.stimAmplitude) > SMALL;

    const bool hasS2 =
        stim.stimPeriodS2 > SMALL
     && stim.nStim2 > 0
     && stim.stimDuration > SMALL
     && mag(stim.stimAmplitude) > SMALL;

    return hasS1 || hasS2;
}

ExternalStimulusProtocol stimulusIO::loadExternalStimulusProtocol
(
    const dictionary& dict
)
{
    ExternalStimulusProtocol stim;

    const dictionary* stimDictPtr = externalStimulusDict(dict);
    if (!stimDictPtr)
    {
        return stim;
    }

    const dictionary& stimDict = *stimDictPtr;
    readMonodomainStimulusBoxes(stimDict, stim.boxes);

    stim.startTimes.setSize(stim.boxes.size());
    if (stimDict.found("stimulusStartTimeList"))
    {
        stimDict.lookup("stimulusStartTimeList") >> stim.startTimes;
        checkSize
        (
            stim.startTimes,
            stim.boxes.size(),
            "stimulusStartTimeList"
        );
    }
    else
    {
        const scalar startTime =
            stimDict.lookupOrDefault<scalar>("stimulusStartTime", 0.0);
        forAll(stim.startTimes, i)
        {
            stim.startTimes[i] = startTime;
        }
    }

    stim.durations = List<scalar>(stim.boxes.size(), 0.0);
    if (stimDict.found("stimulusDurationList"))
    {
        stimDict.lookup("stimulusDurationList") >> stim.durations;
        checkSize
        (
            stim.durations,
            stim.boxes.size(),
            "stimulusDurationList"
        );
    }
    else
    {
        const dimensionedScalar stimulusDuration
        (
            "stimulusDuration", dimTime, stimDict
        );
        forAll(stim.durations, i)
        {
            stim.durations[i] = stimulusDuration.value();
        }
    }
    checkNonNegative(stim.durations, "stimulusDuration");

    stim.intensities = List<scalar>(stim.boxes.size(), 0.0);
    if (stimDict.found("stimulusIntensityList"))
    {
        stimDict.lookup("stimulusIntensityList") >> stim.intensities;
        checkSize
        (
            stim.intensities,
            stim.boxes.size(),
            "stimulusIntensityList"
        );
    }
    else
    {
        const dimensionedScalar stimulusIntensity
        (
            "stimulusIntensity", dimCurrent/dimVolume, stimDict
        );
        forAll(stim.intensities, i)
        {
            stim.intensities[i] = stimulusIntensity.value();
        }
    }

    return stim;
}

word stimulusIO::protocolSuffix(const dictionary& dict)
{
    const dictionary* stimDictPtr = singleCellStimulusDict(dict);
    if (!stimDictPtr)
    {
        return "noStim";
    }
    const dictionary& stimDict = *stimDictPtr;
    const scalar s1 = stimDict.lookupOrDefault<scalar>("stim_period_S1", 0.0);
    const label n1 = stimDict.lookupOrDefault<label>
    (
        "nstim1",
        s1 > SMALL ? 1 : 0
    );
    const label n2 = stimDict.lookupOrDefault<label>("nstim2", 0);

    if (s1 <= SMALL || n1 <= 0)
    {
        return "noStim";
    }

    word out = "S1_" + Foam::name(s1);

    if (n2 > 0)
    {
        const scalar s2 =
            stimDict.lookupOrDefault<scalar>("stim_period_S2", 0.0);
        out += "_S2_" + Foam::name(s2);
    }

    return out;
}
} // namespace Foam
