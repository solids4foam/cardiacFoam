#ifndef ionicModelIO_H
#define ionicModelIO_H

#include "OFstream.H"
#include "PtrList.H"
#include "scalarField.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"

namespace Foam {

class ionicModelIO
{
public:


    // Existing generic writer 
    typedef scalar (*VmTransform)(const scalarField&);

    enum ColumnSource
    {
        COL_VM = 0,
        COL_STATE = 1,
        COL_ALGEBRAIC = 2,
        COL_RATE = 3
    };

    struct SelectionPlan
    {
        wordList names;
        List<label> source;
        List<label> index;
    };

    struct SelectedMapCache
    {
        bool valid = false;
        const void* stateNamesPtr = nullptr;
        const void* algNamesPtr = nullptr;
        label nStates = -1;
        label nAlg = -1;
        wordList exportedNames;
        SelectionPlan plan;
    };

    struct FullPlanCache
    {
        bool valid = false;
        const void* stateNamesPtr = nullptr;
        const void* algNamesPtr = nullptr;
        label nStates = -1;
        label nAlg = -1;
        SelectionPlan plan;
    };

    struct ExportedNamesCache
    {
        bool valid = false;
        const void* stateNamesPtr = nullptr;
        const void* algNamesPtr = nullptr;
        label nStates = -1;
        label nAlg = -1;
        wordList requested;
        wordList filtered;
    };
    
    static void writeHeader
    (
        OFstream& os,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        FullPlanCache& fullPlan
    );
    static void writeSelectedHeader
    (
    OFstream& os,
    const wordList& exportedNames
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

    static void writeSelected
    (
    const scalar t,
    OFstream& os,
    const PtrList<scalarField>& STATES,
    const PtrList<scalarField>& ALGEBRAIC,
    const wordList& exportedNames,
    const char* const stateNames[],
    label nStates,
    const char* const algNames[],
    label nAlg,
    SelectedMapCache& selectedPlan,
    const PtrList<scalarField>& RATES,
    VmTransform transformVm = nullptr
    );

    // Determine whether this timestep should be written (protocol-based)
    static bool shouldWriteStep
    (
        scalar tBegin,
        scalar tEnd,
        const dictionary& dict,
        bool utilitiesMode
    );


    static const wordList& exportedFieldNamesRef
    (
        const wordList& userList,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg,
        ExportedNamesCache& cache
    );

    static void exportStateFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES,
        const wordList& exportedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        SelectedMapCache& selectedPlan,
        PtrList<volScalarField>& outFields
    );
    static void debugPrintFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES,
        const wordList& names,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        SelectedMapCache& selectedPlan,
        label cellI,
        scalar t1,
        scalar t2,
        scalar step
    );


    // Map variable names to Vm/STATE/ALGEBRAIC/RATE indices
    static void mapVariableNames
    (
    const wordList& names,
    const char* const stateNames[],
    int nStates,
    const char* const algNames[],
    int nAlg,
    List<label>& stateIndex,
    List<label>& algIndex,
    List<label>& rateIndex
    );

    static void writeOneSweepRow
    (
        OFstream& os,
        scalar V,
        const wordList& deps,
        const scalarField& STATES,
        const scalarField& ALG,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg,
        const scalarField& RATES,
        SelectedMapCache& selectedPlan
    );

    static void writeSweepHeader
    (
    OFstream& os,
    const wordList& deps
    );

    

    

};

}

#endif
