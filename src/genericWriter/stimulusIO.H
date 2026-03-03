#ifndef stimulusIO_H
#define stimulusIO_H


#include "scalarField.H"
#include "dictionary.H"
#include "word.H"
#include "OFstream.H"
#include "boundBox.H"

namespace Foam
{

class StimulusProtocol
{
public:
    scalar stimStart;
    scalar stimPeriodS1;
    scalar stimDuration;
    scalar stimAmplitude;
    label nStim1;
    scalar stimPeriodS2;
    label nStim2;

    StimulusProtocol()
    :
        stimStart(0.0),
        stimPeriodS1(0.0),
        stimDuration(0.0),
        stimAmplitude(0.0),
        nStim1(0),
        stimPeriodS2(0.0),
        nStim2(0)
    {}
};

class MonodomainStimulusProtocol
{
public:
    List<boundBox> boxes;
    List<scalar> startTimes;
    List<scalar> durations;
    List<scalar> intensities;

    MonodomainStimulusProtocol()
    :
        boxes(0),
        startTimes(0),
        durations(0),
        intensities(0)
    {}
};

class stimulusIO
{
public:

    // Load single-cell S1/S2 protocol from "singleCellStimulus"
    // sub-dictionary. If the sub-dictionary is absent, default to a
    // no-stimulus protocol (all zeros).
    static StimulusProtocol loadStimulusProtocol
    (
        const dictionary& dict
    );

    static void loadStimulusProtocol
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
    );

    static scalar computeStimulus
    (
        scalar t,
        const StimulusProtocol& stim
    );

    static scalar computeStimulus
    (
        scalar t,
        scalar stim_start,
        scalar stim_period_S1,
        scalar stim_duration,
        scalar stim_amplitude,
        label  nstim1,
        scalar stim_period_S2,
        label  nstim2
    );

    // Return true when any S1/S2 pacing branch can generate a non-zero pulse.
    static bool hasActiveStimulus(const StimulusProtocol& stim);

    // Parse monodomain spatial stimulus entries from "monodomainStimulus".
    // If the sub-dictionary is absent, returns an empty protocol.
    static MonodomainStimulusProtocol loadMonodomainStimulusProtocol
    (
        const dictionary& dict
    );

    static word protocolSuffix(const dictionary& dict);

};

} // namespace Foam
#endif
