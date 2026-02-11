#ifndef stimulusIO_H
#define stimulusIO_H


#include "scalarField.H"
#include "dictionary.H"
#include "word.H"
#include "OFstream.H"

namespace Foam
{

class stimulusIO
{
public:

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
        scalar stim_start,
        scalar stim_period_S1,
        scalar stim_duration,
        scalar stim_amplitude,
        label  nstim1,
        scalar stim_period_S2,
        label  nstim2
    );


    // Determine whether we should output this timestep (protocol-based)
    static bool shouldWriteStep
    (
        scalar tBegin,
        scalar tEnd,
        const dictionary& dict,
        bool utilitiesMode
    );


    static word protocolSuffix(const dictionary& dict);

};

} // namespace Foam
#endif
