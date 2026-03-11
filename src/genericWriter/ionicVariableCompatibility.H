#ifndef ionicVariableCompatibility_H
#define ionicVariableCompatibility_H

#include "labelList.H"
#include "wordList.H"

namespace Foam {

class ionicVariableCompatibility
{
public:

    // Used by cache checks in IO planners.
    static bool sameWordList
    (
        const wordList& a,
        const wordList& b
    );

    // Accepts Vm aliases (case/separator-insensitive).
    static bool isVmLikeName(const word& name);

    // Find membrane-voltage state index from model state names.
    // Returns -1 if not found.
    static label findVmStateIndex
    (
        const char* const stateNames[],
        label nStates
    );

    // Find intracellular-calcium state index from model state names.
    // Returns -1 if not found.
    static label findCaiStateIndex
    (
        const char* const stateNames[],
        label nStates
    );

    // Accepts both RATES_<state> and rates_<state>.
    static bool isRateNameRelaxed
    (
        const word& name,
        word& stateName
    );

    // Resolve one requested name against state/algebraic/rate spaces.
    static bool resolveVariable
    (
        const word& name,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg,
        bool& isVm,
        label& stateIdx,
        label& algIdx,
        label& rateIdx
    );

    // Resolve many names into state/algebraic/rate index vectors.
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
};

}

#endif
