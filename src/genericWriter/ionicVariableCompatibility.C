#include "ionicVariableCompatibility.H"
#include "error.H"
#include <string>

namespace Foam
{

namespace
{
    static const char* const ratePrefix = "RATES_";
    static const label ratePrefixLen = 6;

    inline char toLowerAscii(const char c)
    {
        return (c >= 'A' && c <= 'Z') ? static_cast<char>(c - 'A' + 'a') : c;
    }

    std::string canonicalName(const char* raw)
    {
        std::string out;

        for (const char* p = raw; *p; ++p)
        {
            const char c = toLowerAscii(*p);

            // CellML/codegen conventions vary between separators/casing.
            if (c == '_' || c == '-' || c == ' ')
            {
                continue;
            }

            out.push_back(c);
        }

        // AV_<name> and <name> should resolve to the same algebraic symbol.
        if (out.size() > 2 && out[0] == 'a' && out[1] == 'v')
        {
            return out.substr(2);
        }

        return out;
    }

    bool sameCanonicalName(const word& lhs, const char* rhs)
    {
        return canonicalName(lhs.c_str()) == canonicalName(rhs);
    }

    bool isRateNameStrict(const word& name, word& stateName)
    {
        if (name.size() <= ratePrefixLen)
        {
            return false;
        }

        for (label i = 0; i < ratePrefixLen; ++i)
        {
            if (name[i] != ratePrefix[i])
            {
                return false;
            }
        }

        stateName = word(name.c_str() + ratePrefixLen);
        return stateName.size() > 0;
    }
}


bool ionicVariableCompatibility::sameWordList
(
    const wordList& a,
    const wordList& b
)
{
    if (a.size() != b.size())
    {
        return false;
    }

    forAll(a, i)
    {
        if (a[i] != b[i])
        {
            return false;
        }
    }

    return true;
}


bool ionicVariableCompatibility::isVmLikeName(const word& name)
{
    return canonicalName(name.c_str()) == "vm";
}

label ionicVariableCompatibility::findVmStateIndex
(
    const char* const stateNames[],
    const label nStates
)
{
    if (!stateNames || nStates <= 0)
    {
        return -1;
    }

    for (label i = 0; i < nStates; ++i)
    {
        const std::string name = canonicalName(stateNames[i]);
        if
        (
               name == "vm"
            || name == "membranev"
            || name == "voltage"
            || name == "membranepotential"
            || name == "cellv"
        )
        {
            return i;
        }
    }

    // Fallback for models that expose only a single-symbol voltage state.
    if (canonicalName(stateNames[0]) == "v")
    {
        return 0;
    }

    return -1;
}

label ionicVariableCompatibility::findCaiStateIndex
(
    const char* const stateNames[],
    const label nStates
)
{
    if (!stateNames || nStates <= 0)
    {
        return -1;
    }

    for (label i = 0; i < nStates; ++i)
    {
        const std::string name = canonicalName(stateNames[i]);

        if
        (
               name == "cai"
            || name == "calciumcai"
            || name == "intracellularcai"
            || name == "intracellularcalcium"
            || name.find("cai") != std::string::npos
        )
        {
            return i;
        }
    }

    return -1;
}

bool ionicVariableCompatibility::isRateNameRelaxed
(
    const word& name,
    word& stateName
)
{
    if (isRateNameStrict(name, stateName))
    {
        return true;
    }

    if (name.size() <= ratePrefixLen)
    {
        return false;
    }

    for (label i = 0; i < ratePrefixLen; ++i)
    {
        if (toLowerAscii(name[i]) != toLowerAscii(ratePrefix[i]))
        {
            return false;
        }
    }

    stateName = word(name.c_str() + ratePrefixLen);
    return stateName.size() > 0;
}


bool ionicVariableCompatibility::resolveVariable
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
)
{
    isVm = false;
    stateIdx = -1;
    algIdx = -1;
    rateIdx = -1;

    if (name == "Vm" || isVmLikeName(name))
    {
        isVm = true;
        stateIdx = 0;
        return true;
    }

    // Match STATE variable
    for (label s = 0; s < nStates; ++s)
    {
        if (name == stateNames[s] || sameCanonicalName(name, stateNames[s]))
        {
            stateIdx = s;
            return true;
        }
    }

    // Match ALGEBRAIC variable
    for (label a = 0; a < nAlg; ++a)
    {
        if (name == algNames[a] || sameCanonicalName(name, algNames[a]))
        {
            algIdx = a;
            return true;
        }
    }

    // Match RATE variable using RATES_<stateName>
    {
        word rateStateName;
        if (isRateNameRelaxed(name, rateStateName))
        {
            for (label s = 0; s < nStates; ++s)
            {
                if
                (
                    rateStateName == stateNames[s]
                 || sameCanonicalName(rateStateName, stateNames[s])
                )
                {
                    rateIdx = s;
                    return true;
                }
            }
        }
    }

    return false;
}


void ionicVariableCompatibility::mapVariableNames
(
    const wordList& names,
    const char* const stateNames[],
    int nStates,
    const char* const algNames[],
    int nAlg,
    List<label>& stateIndex,
    List<label>& algIndex,
    List<label>& rateIndex
)
{
    stateIndex.setSize(names.size(), -1);
    algIndex.setSize(names.size(), -1);
    rateIndex.setSize(names.size(), -1);

    forAll(names, k)
    {
        bool isVm = false;
        label sIdx = -1;
        label aIdx = -1;
        label rIdx = -1;

        if
        (
            resolveVariable
            (
                names[k],
                stateNames,
                nStates,
                algNames,
                nAlg,
                isVm,
                sIdx,
                aIdx,
                rIdx
            )
        )
        {
            if (sIdx >= 0)
            {
                stateIndex[k] = sIdx;
            }
            if (aIdx >= 0)
            {
                algIndex[k] = aIdx;
            }
            if (rIdx >= 0)
            {
                rateIndex[k] = rIdx;
            }
            continue;
        }

        FatalErrorInFunction
            << "Unknown debug variable: " << names[k]
            << ". Expected Vm, one of the STATE names, one of the ALGEBRAIC "
            << "names, or a rate name of the form RATES_<stateName>."
            << exit(FatalError);
    }
}

} // namespace Foam
