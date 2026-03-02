#ifndef stimulusIO_H
#define stimulusIO_H


#include "scalarField.H"
#include "dictionary.H"
#include "word.H"
#include "OFstream.H"

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

class stimulusIO
{
public:

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
    

    static word protocolSuffix(const dictionary& dict);

};

} // namespace Foam
#endif
