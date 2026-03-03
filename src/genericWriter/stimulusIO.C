#include "stimulusIO.H"
#include "IOstreams.H"
#include <cmath>

namespace Foam
{
StimulusProtocol stimulusIO::loadStimulusProtocol
(
    const dictionary& dict
)
{
    StimulusProtocol stim;

    // Required keys
    const char* required[] =
    {
        "stim_start",
        "stim_period_S1",
        "stim_duration",
        "stim_amplitude"
    };

    for (const char* k : required)
    {
        if (!dict.found(k))
        {
            FatalErrorInFunction
                << "Missing required stimulus key: " << k
                << exit(FatalError);
        }
    }

    stim.stimStart = readScalar(dict.lookup("stim_start"));
    stim.stimPeriodS1 = readScalar(dict.lookup("stim_period_S1"));
    stim.stimDuration = readScalar(dict.lookup("stim_duration"));
    stim.stimAmplitude = readScalar(dict.lookup("stim_amplitude"));

    // Defaults for optional keys
    stim.nStim1 = 1;
    stim.nStim2 = 0;
    stim.stimPeriodS2 = 0.0;

    if (dict.found("nstim1"))
    {
        stim.nStim1 = readLabel(dict.lookup("nstim1"));
    }

    if (dict.found("stim_period_S2"))
    {
        stim.stimPeriodS2 = readScalar(dict.lookup("stim_period_S2"));
    }

    if (dict.found("nstim2"))
    {
        stim.nStim2 = readLabel(dict.lookup("nstim2"));
    }

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

    word stimulusIO::protocolSuffix(const dictionary& dict)
    {
        scalar s1 = readScalar(dict.lookup("stim_period_S1"));
        scalar n2 = dict.found("nstim2") ? readScalar(dict.lookup("nstim2")) : 0;

        word out = "S1_" + Foam::name(s1);

        if (n2 > 0)
        {
            scalar s2 = readScalar(dict.lookup("stim_period_S2"));
            out += "_S2_" + Foam::name(s2);
        }

        return out;
    }
} // namespace Foam

