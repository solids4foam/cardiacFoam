#ifndef ionicModelIO_H
#define ionicModelIO_H

#include "OFstream.H"
#include "PtrList.H"
#include "scalarField.H"
#include "volFields.H"

namespace Foam {

class ionicModelIO
{
public:

    static void loadStimulusConstants
    (
        const dictionary& dict,
        scalarField& CONSTANTS,

        // model index positions
        label stim_start,
        label stim_period_S1,
        label stim_duration,
        label stim_amplitude,
        label nstim1,
        label stim_period_S2,
        label nstim2
    );


    // Function pointer for model-specific Vm conversion
    typedef scalar (*VmTransform)(const scalarField&);

    // Write header line
    static void writeHeader
    (
        OFstream& os,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg
    );

    // Generic writer method
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

};
}

#endif
