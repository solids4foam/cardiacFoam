#ifndef ionicModelIO_H
#define ionicModelIO_H

#include "OFstream.H"
#include "PtrList.H"
#include "scalarField.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "ionicModel.H"

namespace Foam {

class ionicModelIO
{
public:

    // Create output filename based on model/tissue/S1-S2 protocol
    static fileName createOutputFile
    (
        const dictionary& dict,
        const word& modelName,
        const word& tissueName,
        const Time& runTime
    );

    // Decide if this timestep should be written
    static bool shouldWriteStep
    (
        const dictionary& dict,
        const scalar currentTime
    );

    // Write a timestep (handles protocol logic)
    static void writeTimestep
    (
        ionicModel& model,
        OFstream& os,
        const dictionary& dict,
        const Time& runTime
    );

    // load stimulus constants
    static void loadStimulusConstants
    (
        const dictionary& dict,
        scalarField& CONSTANTS,

        label stim_start,
        label stim_period_S1,
        label stim_duration,
        label stim_amplitude,
        label nstim1,
        label stim_period_S2,
        label nstim2
    );

    // Existing generic writer
    typedef scalar (*VmTransform)(const scalarField&);

    static void writeHeader
    (
        OFstream& os,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg
    );
    static void write
    (
        scalar t,
        OFstream& os,
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES,
        VmTransform transformVm
    );

    static wordList exportedFieldNames
    (
        const wordList& userList,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg
    );

    static void exportStateFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const wordList& exportedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        List<volScalarField*>& outFields
    );
    static void debugPrintFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const wordList& names,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        label cellI,
        scalar t1,
        scalar t2,
        scalar step
    );


    // Map variable names to STATE or ALGEBRAIC indices
    static void mapVariableNames
    (
    const wordList& names,
    const char* const stateNames[],
    int nStates,
    const char* const algNames[],
    int nAlg,
    List<label>& stateIndex,
    List<label>& algIndex
    );



};

}

#endif


