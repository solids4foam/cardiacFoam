#include "stimulusIO.H"
#include "IOstreams.H"
#include <cmath>

namespace Foam
{

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

    C[stim_start]     = readScalar(dict.lookup("stim_start"));
    C[stim_period_S1] = readScalar(dict.lookup("stim_period_S1"));
    C[stim_duration]  = readScalar(dict.lookup("stim_duration"));
    C[stim_amplitude] = readScalar(dict.lookup("stim_amplitude"));

    // Defaults for optional keys
    C[nstim1] = 1;          // default: one S1 pulse train period count
    C[nstim2] = 0;          // default: no S2 unless requested
    C[stim_period_S2] = 0;  // default: disabled

if (dict.found("nstim1"))
    C[nstim1] = readScalar(dict.lookup("nstim1"));

if (dict.found("stim_period_S2"))
    C[stim_period_S2] = readScalar(dict.lookup("stim_period_S2"));

if (dict.found("nstim2"))
    C[nstim2] = readScalar(dict.lookup("nstim2"));

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







    bool stimulusIO::shouldWriteStep
    (
        scalar tBegin,
        scalar tEnd,
        const dictionary& dict,
        bool utilitiesMode
    )
    {
        // For utilities (like sweepCurrents): always write
        if (utilitiesMode)
            return true;

        // User controls when to start writing
        scalar writeAfterTime = 0.0;
        if (dict.found("writeAfterTime"))
        {
            writeAfterTime = readScalar(dict.lookup("writeAfterTime"));
        }
        // Start writing once we pass that time
        return (tEnd >= writeAfterTime);
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



