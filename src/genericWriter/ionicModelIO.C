/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ionicModelIO.H"
#include "ionicVariableCompatibility.H"
#include "DynamicList.H"

#include <string>

namespace Foam
{

    namespace
    {
        inline char toLowerAscii(const char c)
        {
            return (c >= 'A' && c <= 'Z') ? static_cast<char>(c - 'A' + 'a') : c;
        }

        word suffixScalar(const scalar value)
        {
            const std::string raw(Foam::name(value).c_str());
            std::string out;
            out.reserve(raw.size());

            for (const char c : raw)
            {
                if (c == '.')
                {
                    out.push_back('p');
                }
                else if (c == '-')
                {
                    out.push_back('m');
                }
                else if (c == '+')
                {
                    out.push_back('p');
                }
                else
                {
                    out.push_back(c);
                }
            }

            return word(out);
        }

        std::string canonicalName(const char* raw)
        {
            std::string out;

            for (const char* p = raw; *p; ++p)
            {
                const char c = toLowerAscii(*p);

                if (c == '_' || c == '-' || c == ' ')
                {
                    continue;
                }

                out.push_back(c);
            }

            if (out.size() > 2 && out[0] == 'a' && out[1] == 'v')
            {
                return out.substr(2);
            }

            return out;
        }

        label findConstantIndex
        (
            const word& requestedName,
            const char* const constantNames[],
            const label nConstants
        )
        {
            if (!constantNames || nConstants <= 0)
            {
                return -1;
            }

            for (label i = 0; i < nConstants; ++i)
            {
                if (requestedName == constantNames[i])
                {
                    return i;
                }
            }

            const std::string canonicalRequested =
                canonicalName(requestedName.c_str());

            for (label i = 0; i < nConstants; ++i)
            {
                if (canonicalRequested == canonicalName(constantNames[i]))
                {
                    return i;
                }
            }

            return -1;
        }

        wordList constantNamesList
        (
            const char* const constantNames[],
            const label nConstants
        )
        {
            wordList names(max(nConstants, label(0)));
            for (label i = 0; i < nConstants; ++i)
            {
                names[i] = word(constantNames[i]);
            }
            return names;
        }

        void validateConstantMetadata
        (
            const char* const constantNames[],
            const label nConstants,
            const word& modelName
        )
        {
            if (!constantNames || nConstants <= 0)
            {
                FatalErrorInFunction
                    << "Constant metadata is not available for ionic model "
                    << modelName << "." << exit(FatalError);
            }

            for (label i = 0; i < nConstants; ++i)
            {
                if (!constantNames[i] || constantNames[i][0] == '\0')
                {
                    FatalErrorInFunction
                        << "Invalid constant metadata for ionic model "
                        << modelName << ": constant index " << i
                        << " has no name. Fix the model's CONSTANTS_NAMES "
                        << "array so every NUM_CONSTANTS entry has a real "
                        << "semantic name."
                        << exit(FatalError);
                }
            }
        }

        struct ConstantOverrideOp
        {
            DynamicList<label> indices;
            DynamicList<scalar> values;
            DynamicList<word> requestedNames;
        };

        word tissueOverrideScopeName(const label tissueFlag)
        {
            return tissueFlag == 1 ? word("epicardialCells")
                 : tissueFlag == 2 ? word("mCells")
                 : tissueFlag == 3 ? word("endocardialCells")
                 : tissueFlag == 4 ? word("myocyte")
                                   : word();
        }

        bool isConstantOverrideScopeName(const word& name)
        {
            return name == "global"
                || name == "epicardialCells"
                || name == "mCells"
                || name == "endocardialCells"
                || name == "myocyte";
        }

        bool isKnownOverrideScopeName
        (
            const word& name,
            const wordList& knownScopeNames
        )
        {
            if (isConstantOverrideScopeName(name))
            {
                return true;
            }

            forAll(knownScopeNames, i)
            {
                if (knownScopeNames[i] == name)
                {
                    return true;
                }
            }

            return false;
        }

        void collectConstantOverrideOp
        (
            const dictionary& opDict,
            const word& scopeName,
            const word& opName,
            const char* const constantNames[],
            const label nConstants,
            const word& modelName,
            List<label>& seen,
            const List<label>& otherSeen,
            ConstantOverrideOp& op
        )
        {
            forAllConstIter(dictionary, opDict, iter)
            {
                const entry& e = iter();
                const word requestedName(e.keyword());

                if (e.isDict())
                {
                    FatalErrorInFunction
                        << "ionicConstantOverrides." << scopeName
                        << "." << opName
                        << " entry '" << requestedName
                        << "' for ionic model " << modelName
                        << " must be a scalar value, not a dictionary."
                        << exit(FatalError);
                }

                const label constantI =
                    findConstantIndex(requestedName, constantNames, nConstants);

                if (constantI < 0)
                {
                    FatalErrorInFunction
                        << "Unknown ionic constant '" << requestedName
                        << "' in ionicConstantOverrides." << scopeName
                        << "." << opName
                        << " for ionic model " << modelName << "." << nl
                        << "Available constants: "
                        << constantNamesList(constantNames, nConstants)
                        << exit(FatalError);
                }

                if (seen[constantI])
                {
                    FatalErrorInFunction
                        << "Duplicate ionic constant override for '"
                        << constantNames[constantI]
                        << "' in ionicConstantOverrides." << scopeName
                        << "." << opName
                        << " for ionic model " << modelName << "."
                        << exit(FatalError);
                }

                if (otherSeen[constantI])
                {
                    FatalErrorInFunction
                        << "Ambiguous ionic constant override for '"
                        << constantNames[constantI]
                        << "' in ionic model " << modelName
                        << ": the same constant appears in both scale and set."
                        << exit(FatalError);
                }

                seen[constantI] = 1;
                op.indices.append(constantI);
                op.values.append(readScalar(opDict.lookup(requestedName)));
                op.requestedNames.append(requestedName);
            }
        }

        void applyConstantOverrideScope
        (
            scalarField& constants,
            const char* const constantNames[],
            const label nConstants,
            const word& modelName,
            const dictionary& scopeDict,
            const word& scopeName
        )
        {
            forAllConstIter(dictionary, scopeDict, iter)
            {
                const word entryName(iter().keyword());
                if (entryName != "scale" && entryName != "set")
                {
                    FatalErrorInFunction
                        << "Unsupported ionicConstantOverrides."
                        << scopeName << " entry '" << entryName
                        << "' for ionic model " << modelName << "."
                        << nl
                        << "Supported entries are 'scale' and 'set'."
                        << exit(FatalError);
                }
            }

            List<label> scaleSeen(nConstants, 0);
            List<label> setSeen(nConstants, 0);
            ConstantOverrideOp scaleOp;
            ConstantOverrideOp setOp;

            if (scopeDict.found("scale"))
            {
                collectConstantOverrideOp
                (
                    scopeDict.subDict("scale"),
                    scopeName,
                    "scale",
                    constantNames,
                    nConstants,
                    modelName,
                    scaleSeen,
                    setSeen,
                    scaleOp
                );
            }

            if (scopeDict.found("set"))
            {
                collectConstantOverrideOp
                (
                    scopeDict.subDict("set"),
                    scopeName,
                    "set",
                    constantNames,
                    nConstants,
                    modelName,
                    setSeen,
                    scaleSeen,
                    setOp
                );
            }

            forAll(scaleOp.indices, opI)
            {
                const label constantI = scaleOp.indices[opI];
                const scalar oldValue = constants[constantI];
                constants[constantI] *= scaleOp.values[opI];

                Info<< "ionicConstantOverrides." << scopeName << ": "
                    << modelName << " scaled " << constantNames[constantI]
                    << " (" << scaleOp.requestedNames[opI] << ")"
                    << " by " << scaleOp.values[opI]
                    << ": " << oldValue << " -> " << constants[constantI]
                    << nl;
            }

            forAll(setOp.indices, opI)
            {
                const label constantI = setOp.indices[opI];
                const scalar oldValue = constants[constantI];
                constants[constantI] = setOp.values[opI];

                Info<< "ionicConstantOverrides." << scopeName << ": "
                    << modelName << " set " << constantNames[constantI]
                    << " (" << setOp.requestedNames[opI] << ")"
                    << ": " << oldValue << " -> " << constants[constantI]
                    << nl;
            }
        }

        void emitHeader
        (
            OFstream& os,
            const wordList& names
        )
        {
            os << "time";
            forAll(names, k)
            {
                os << " " << names[k];
            }
            os << nl;
        }

        void emitRow
        (
            const scalar t,
            OFstream& os,
            const scalarField& S,
            const scalarField& A,
            const scalarField& R,
            const ionicModelIO::SelectionPlan& plan,
            Foam::ionicModelIO::VmTransform transformVm
        )
        {
            os << t;

            forAll(plan.source, k)
            {
                const label src = plan.source[k];
                const label idx = plan.index[k];

                if (src == ionicModelIO::COL_VM)
                {
                    os << " " << (transformVm ? transformVm(S) : S[0]);
                }
                else if (src == ionicModelIO::COL_STATE)
                {
                    os << " " << S[idx];
                }
                else if (src == ionicModelIO::COL_ALGEBRAIC)
                {
                    os << " " << A[idx];
                }
                else if (src == ionicModelIO::COL_RATE)
                {
                    os << " " << R[idx];
                }
                else
                {
                    FatalErrorInFunction
                        << "Unexpected selection source id " << src
                        << " in write plan." << exit(FatalError);
                }
            }

            os << nl;
        }

        const ionicModelIO::SelectionPlan& selectedPlan
        (
            const wordList& exportedNames,
            const char* const stateNames[],
            const label nStates,
            const char* const algNames[],
            const label nAlg,
            ionicModelIO::SelectedMapCache& cache
        )
        {
            const bool cacheHit =
                   cache.valid
                && cache.stateNamesPtr == static_cast<const void*>(stateNames)
                && cache.algNamesPtr == static_cast<const void*>(algNames)
                && cache.nStates == nStates
                && cache.nAlg == nAlg
                && ionicVariableCompatibility::sameWordList
                   (
                       cache.exportedNames,
                       exportedNames
                   );

            if (!cacheHit)
            {
                cache.stateNamesPtr = static_cast<const void*>(stateNames);
                cache.algNamesPtr = static_cast<const void*>(algNames);
                cache.nStates = nStates;
                cache.nAlg = nAlg;
                cache.exportedNames = exportedNames;

                List<label> stateIndex;
                List<label> algIndex;
                List<label> rateIndex;
                Foam::ionicModelIO::mapVariableNames
                (
                    exportedNames,
                    stateNames,
                    nStates,
                    algNames,
                    nAlg,
                    stateIndex,
                    algIndex,
                    rateIndex
                );

                cache.plan.names = exportedNames;
                cache.plan.source.setSize(exportedNames.size());
                cache.plan.index.setSize(exportedNames.size());

                forAll(exportedNames, k)
                {
                    if
                    (
                        exportedNames[k] == "Vm"
                     || ionicVariableCompatibility::isVmLikeName(exportedNames[k])
                    )
                    {
                        cache.plan.source[k] = ionicModelIO::COL_VM;
                        cache.plan.index[k] = 0;
                    }
                    else if (stateIndex[k] >= 0)
                    {
                        cache.plan.source[k] = ionicModelIO::COL_STATE;
                        cache.plan.index[k] = stateIndex[k];
                    }
                    else if (rateIndex[k] >= 0)
                    {
                        cache.plan.source[k] = ionicModelIO::COL_RATE;
                        cache.plan.index[k] = rateIndex[k];
                    }
                    else
                    {
                        cache.plan.source[k] = ionicModelIO::COL_ALGEBRAIC;
                        cache.plan.index[k] = algIndex[k];
                    }
                }
                cache.valid = true;
            }

            return cache.plan;
        }

        const ionicModelIO::SelectionPlan& fullPlan
        (
            const char* const stateNames[],
            const label nStates,
            const char* const algNames[],
            const label nAlg,
            ionicModelIO::FullPlanCache& cache
        )
        {
            const bool cacheHit =
                   cache.valid
                && cache.stateNamesPtr == static_cast<const void*>(stateNames)
                && cache.algNamesPtr == static_cast<const void*>(algNames)
                && cache.nStates == nStates
                && cache.nAlg == nAlg;

            if (!cacheHit)
            {
                cache.stateNamesPtr = static_cast<const void*>(stateNames);
                cache.algNamesPtr = static_cast<const void*>(algNames);
                cache.nStates = nStates;
                cache.nAlg = nAlg;

                ionicModelIO::SelectionPlan& plan = cache.plan;
                const label nCols = 1 + (nStates - 1) + nAlg + nStates;
                plan.names.setSize(nCols);
                plan.source.setSize(nCols);
                plan.index.setSize(nCols);

                label k = 0;

                plan.names[k] = "Vm";
                plan.source[k] = ionicModelIO::COL_VM;
                plan.index[k] = 0;
                ++k;

                for (label s = 1; s < nStates; ++s)
                {
                    plan.names[k] = stateNames[s];
                    plan.source[k] = ionicModelIO::COL_STATE;
                    plan.index[k] = s;
                    ++k;
                }

                for (label a = 0; a < nAlg; ++a)
                {
                    plan.names[k] = algNames[a];
                    plan.source[k] = ionicModelIO::COL_ALGEBRAIC;
                    plan.index[k] = a;
                    ++k;
                }

                for (label s = 0; s < nStates; ++s)
                {
                    plan.names[k] = word("RATES_") + word(stateNames[s]);
                    plan.source[k] = ionicModelIO::COL_RATE;
                    plan.index[k] = s;
                    ++k;
                }

                cache.valid = true;
            }

            return cache.plan;
        }
    }

    void Foam::ionicModelIO::writeHeader
    (
        OFstream& os,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        FullPlanCache& fullPlanCache
    )
    {
        emitHeader
        (
            os,
            fullPlan(stateNames, nStates, algNames, nAlg, fullPlanCache).names
        );
    }

    void Foam::ionicModelIO::writeSelectedHeader
    (
        OFstream& os,
        const wordList& exportedNames
    )
    {
        emitHeader(os, exportedNames);
    }


    void Foam::ionicModelIO::write
    (
        scalar t,
        OFstream& os,
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES,
        VmTransform transformVm
    )
    {
        const scalarField& S = STATES[0];
        const scalarField& A = ALGEBRAIC[0];
        const scalarField& R = RATES[0];

        // Keep full-row emission tight; selected/export paths share one plan engine.
        os << t << " " << (transformVm ? transformVm(S) : S[0]);

        for (label i = 1; i < S.size(); ++i)
        {
            os << " " << S[i];
        }

        forAll(A, i)
        {
            os << " " << A[i];
        }

        forAll(R, i)
        {
            os << " " << R[i];
        }

        os << nl;
    }

    void Foam::ionicModelIO::writeSelected
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
        SelectedMapCache& selectedPlanCache,
        const PtrList<scalarField>& RATES,
        VmTransform transformVm
    )
    {
        // 1. mapping (cached between timesteps)
        const SelectionPlan& plan =
            selectedPlan
            (
                exportedNames,
                stateNames,
                nStates,
                algNames,
                nAlg,
                selectedPlanCache
            );

        // 2. single cell only
        const scalarField& S = STATES[0];
        const scalarField& A = ALGEBRAIC[0];
        const scalarField& R = RATES[0];

        emitRow(t, os, S, A, R, plan, transformVm);
    }

    bool Foam::ionicModelIO::shouldWriteStep
    (
        scalar tBegin,
        scalar tEnd,
        const dictionary& dict,
        bool utilitiesMode
    )
    {
        // For utilities (like sweepCurrents): always write
        if (utilitiesMode)
        {
            return true;
        }

        // User controls when to start writing
        scalar writeAfterTime = 0.0;
        if (dict.found("writeAfterTime"))
        {
            writeAfterTime = readScalar(dict.lookup("writeAfterTime"));
        }

        // Start writing once we pass that time
        if (tEnd < writeAfterTime)
        {
            return false;
        }

        scalar writeFrequency = 0.0;
        if (dict.found("writeFrequency"))
        {
            writeFrequency = readScalar(dict.lookup("writeFrequency"));
        }

        if (writeFrequency > 0.0)
        {
            // Write only if we cross a frequency boundary between tBegin and tEnd
            label stepBegin = std::floor(tBegin / writeFrequency);
            label stepEnd = std::floor(tEnd / writeFrequency);
            return (stepEnd > stepBegin);
        }

        return true;
    }

    void Foam::ionicModelIO::applyConstantOverrides
    (
        scalarField& constants,
        const char* const constantNames[],
        const label nConstants,
        const dictionary& dict,
        const word& modelName
    )
    {
        applyConstantOverrides
        (
            constants,
            constantNames,
            nConstants,
            dict,
            modelName,
            label(-1)
        );
    }


    void Foam::ionicModelIO::applyConstantOverrides
    (
        scalarField& constants,
        const char* const constantNames[],
        const label nConstants,
        const dictionary& dict,
        const word& modelName,
        const label tissueFlag
    )
    {
        if (!dict.found("ionicConstantOverrides"))
        {
            return;
        }

        if (constants.empty())
        {
            FatalErrorInFunction
                << "ionicConstantOverrides was requested for ionic model "
                << modelName
                << ", but constant storage is not available."
                << exit(FatalError);
        }

        validateConstantMetadata(constantNames, nConstants, modelName);

        if (constants.size() != nConstants)
        {
            FatalErrorInFunction
                << "ionicConstantOverrides for ionic model " << modelName
                << " found " << nConstants << " constant names but "
                << constants.size() << " stored constant values."
                << exit(FatalError);
        }

        const dictionary& overrides = dict.subDict("ionicConstantOverrides");

        forAllConstIter(dictionary, overrides, iter)
        {
            const word entryName(iter().keyword());
            if (!isConstantOverrideScopeName(entryName))
            {
                // Not one of the fixed global/anatomical/myocyte scopes this
                // tissueFlag-based overload understands. This may be a
                // namedRegions-only scope (e.g. a scar/disease region name),
                // which is validated and applied separately by the
                // word-scoped applyConstantOverrides overload used by
                // ionicModel::constantsForRegion(). Skip it here rather
                // than fatal, since this overload has no visibility into
                // which named regions are legitimately declared elsewhere
                // in the dictionary.
                continue;
            }

            if (!iter().isDict())
            {
                FatalErrorInFunction
                    << "ionicConstantOverrides entry '" << entryName
                    << "' for ionic model " << modelName
                    << " must be a dictionary."
                    << exit(FatalError);
            }
        }

        if (overrides.found("global"))
        {
            applyConstantOverrideScope
            (
                constants,
                constantNames,
                nConstants,
                modelName,
                overrides.subDict("global"),
                "global"
            );
        }

        const word tissueScopeName = tissueOverrideScopeName(tissueFlag);
        if (!tissueScopeName.empty() && overrides.found(tissueScopeName))
        {
            applyConstantOverrideScope
            (
                constants,
                constantNames,
                nConstants,
                modelName,
                overrides.subDict(tissueScopeName),
                tissueScopeName
            );
        }
    }


    void Foam::ionicModelIO::applyConstantOverrides
    (
        scalarField& constants,
        const char* const constantNames[],
        const label nConstants,
        const dictionary& dict,
        const word& modelName,
        const word& scopeName,
        const word& baselineScopeName,
        const wordList& knownScopeNames
    )
    {
        if (!dict.found("ionicConstantOverrides"))
        {
            return;
        }

        if (constants.empty())
        {
            FatalErrorInFunction
                << "ionicConstantOverrides was requested for ionic model "
                << modelName
                << ", but constant storage is not available."
                << exit(FatalError);
        }

        validateConstantMetadata(constantNames, nConstants, modelName);

        if (constants.size() != nConstants)
        {
            FatalErrorInFunction
                << "ionicConstantOverrides for ionic model " << modelName
                << " found " << nConstants << " constant names but "
                << constants.size() << " stored constant values."
                << exit(FatalError);
        }

        const dictionary& overrides = dict.subDict("ionicConstantOverrides");

        forAllConstIter(dictionary, overrides, iter)
        {
            const word entryName(iter().keyword());
            if (!isKnownOverrideScopeName(entryName, knownScopeNames))
            {
                FatalErrorInFunction
                    << "Unsupported ionicConstantOverrides entry '"
                    << entryName << "' for ionic model " << modelName << "."
                    << nl
                    << "Declared regions: " << knownScopeNames
                    << " (plus global/anatomical/myocyte scopes)."
                    << exit(FatalError);
            }

            if (!iter().isDict())
            {
                FatalErrorInFunction
                    << "ionicConstantOverrides entry '" << entryName
                    << "' for ionic model " << modelName
                    << " must be a dictionary."
                    << exit(FatalError);
            }
        }

        // constantsForRegion() starts from constantsForTissue(), which has
        // already applied global and the explicit baseline scope. Layer the
        // region's own scope only when it names a distinct scope; validation
        // above still catches typos in either fixed or declared scopes.
        if
        (
            !scopeName.empty()
         && scopeName != baselineScopeName
         && overrides.found(scopeName)
        )
        {
            applyConstantOverrideScope
            (
                constants, constantNames, nConstants, modelName,
                overrides.subDict(scopeName), scopeName
            );
        }
    }


    Foam::word Foam::ionicModelIO::constantOverrideOutputSuffix
    (
        const dictionary& dict
    )
    {
        if (dict.found("outputSuffix"))
        {
            return dict.lookupOrDefault<word>("outputSuffix", word());
        }

        if (!dict.found("ionicConstantOverrides"))
        {
            return word();
        }

        const dictionary& overrides = dict.subDict("ionicConstantOverrides");
        word suffix;
        wordList scopeNames(5);
        scopeNames[0] = "global";
        scopeNames[1] = "endocardialCells";
        scopeNames[2] = "mCells";
        scopeNames[3] = "epicardialCells";
        scopeNames[4] = "myocyte";

        wordList opNames(2);
        opNames[0] = "scale";
        opNames[1] = "set";

        forAll(scopeNames, scopeI)
        {
            const word& scopeName = scopeNames[scopeI];
            if (!overrides.found(scopeName))
            {
                continue;
            }

            const dictionary& scopeDict = overrides.subDict(scopeName);
            forAll(opNames, opI)
            {
                const word& opName = opNames[opI];

                if (!scopeDict.found(opName))
                {
                    continue;
                }

                const dictionary& opDict = scopeDict.subDict(opName);

                forAllConstIter(dictionary, opDict, iter)
                {
                    const entry& e = iter();
                    if (e.isDict())
                    {
                        continue;
                    }

                    if (!suffix.empty())
                    {
                        suffix += "_";
                    }

                    const word constantName(e.keyword());
                    if (scopeName != "global")
                    {
                        suffix += scopeName;
                        suffix += "_";
                    }
                    suffix += opName;
                    suffix += "_";
                    suffix += constantName;
                    suffix += "_";
                    suffix +=
                        suffixScalar(readScalar(opDict.lookup(constantName)));
                }
            }
        }

        return suffix;
    }


    const Foam::wordList& Foam::ionicModelIO::exportedFieldNamesRef
    (
        const wordList& userList,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg,
        ExportedNamesCache& cache
    )
    {
        if (!userList.size())
        {
            static const wordList empty;
            return empty;
        }

        const bool cacheHit =
               cache.valid
            && cache.stateNamesPtr == static_cast<const void*>(stateNames)
            && cache.algNamesPtr == static_cast<const void*>(algNames)
            && cache.nStates == nStates
            && cache.nAlg == nAlg
            && ionicVariableCompatibility::sameWordList(cache.requested, userList);
        if (cacheHit)
        {
            return cache.filtered;
        }

        wordList filtered;
        forAll(userList, k)
        {
            const word& name = userList[k];
            bool isVm = false;
            label stateIdx = -1;
            label algIdx = -1;
            label rateIdx = -1;
            const bool found =
                ionicVariableCompatibility::resolveVariable
                (
                    name,
                    stateNames,
                    nStates,
                    algNames,
                    nAlg,
                    isVm,
                    stateIdx,
                    algIdx,
                    rateIdx
                )
             && (isVm || stateIdx >= 0 || algIdx >= 0 || rateIdx >= 0);

            if (found)
            {
                filtered.append(name);
            }
            else
            {
                WarningInFunction
                    << "Ignoring unknown ionic model variable '" << name
                    << "'." << nl;
            }
        }

        cache.valid = true;
        cache.stateNamesPtr = static_cast<const void*>(stateNames);
        cache.algNamesPtr = static_cast<const void*>(algNames);
        cache.nStates = nStates;
        cache.nAlg = nAlg;
        cache.requested = userList;
        cache.filtered = filtered;

        return cache.filtered;
    }

    void Foam::ionicModelIO::exportStateFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES,
        const wordList& exportedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        SelectedMapCache& selectedPlanCache,
        PtrList<volScalarField>& outFields
    )
    {
        const SelectionPlan& plan =
            selectedPlan
            (
                exportedNames,
                stateNames,
                nStates,
                algNames,
                nAlg,
                selectedPlanCache
            );

        if (outFields.size() != plan.source.size())
        {
            FatalErrorInFunction
                << "Mismatch between selected export variables ("
                << plan.source.size() << ") and allocated output fields ("
                << outFields.size() << ")."
                << exit(FatalError);
        }

        // 2. Populate volScalarFields
        forAll(STATES, cellI)
        {
            const scalarField& S = STATES[cellI];
            const scalarField& A = ALGEBRAIC[cellI];
            const scalarField& R = RATES[cellI];

            forAll(outFields, k)
            {
                const label src = plan.source[k];
                const label idx = plan.index[k];

                if (src == ionicModelIO::COL_VM || src == ionicModelIO::COL_STATE)
                {
                    outFields[k][cellI] = S[idx];
                }
                else if (src == ionicModelIO::COL_RATE)
                {
                    outFields[k][cellI] = R[idx];
                }
                else if (src == ionicModelIO::COL_ALGEBRAIC)
                {
                    outFields[k][cellI] = A[idx];
                }
                else
                {
                    FatalErrorInFunction
                        << "Unexpected selection source id " << src
                        << " in export plan." << exit(FatalError);
                }
            }
        }

        // 3. Boundaries
        forAll(outFields, k)
        {
            outFields[k].correctBoundaryConditions();
        }
    }

    void Foam::ionicModelIO::importStateFields
    (
        PtrList<scalarField>& STATES,
        const volScalarField& Vm,
        const PtrList<volScalarField>& inFields,
        const wordList& fieldNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        SelectedMapCache& selectedPlanCache
    )
    {
        const SelectionPlan& plan =
            selectedPlan
            (
                fieldNames,
                stateNames,
                nStates,
                algNames,
                nAlg,
                selectedPlanCache
            );

        if (inFields.size() != plan.source.size())
        {
            FatalErrorInFunction
                << "Mismatch between selected import variables ("
                << plan.source.size() << ") and provided input fields ("
                << inFields.size() << ")."
                << exit(FatalError);
        }

        if (STATES.size() != Vm.size())
        {
            FatalErrorInFunction
                << "Vm field size " << Vm.size()
                << " does not match ionic state count " << STATES.size() << "."
                << exit(FatalError);
        }

        const scalarField& VmValues = Vm.primitiveField();
        const label vmStateIndex =
            ionicVariableCompatibility::findVmStateIndex(stateNames, nStates);

        forAll(STATES, cellI)
        {
            scalarField& S = STATES[cellI];

            // Keep the ionic-model voltage state synchronized with the
            // externally provided Vm field even when the caller only imports
            // auxiliary manufactured states such as u1/u2/u3.
            if (vmStateIndex >= 0)
            {
                S[vmStateIndex] = VmValues[cellI];
            }

            forAll(inFields, k)
            {
                const label src = plan.source[k];
                const label idx = plan.index[k];

                if (src == ionicModelIO::COL_VM)
                {
                    S[idx] = VmValues[cellI];
                }
                else if (src == ionicModelIO::COL_STATE)
                {
                    S[idx] = inFields[k][cellI];
                }
                else
                {
                    FatalErrorInFunction
                        << "Cannot import field '" << fieldNames[k]
                        << "' into ionic states because it resolves to source "
                        << src
                        << ". Only Vm/state variables are supported for import."
                        << exit(FatalError);
                }
            }
        }
    }

    void Foam::ionicModelIO::debugPrintFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES,
        const wordList& printedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        SelectedMapCache& selectedPlanCache,
        label cellI,
        scalar t1,
        scalar t2,
        scalar step
    )
    {
        if (printedNames.empty())
            return;

        const SelectionPlan& plan =
            selectedPlan
            (
                printedNames,
                stateNames,
                nStates,
                algNames,
                nAlg,
                selectedPlanCache
            );

        const scalarField& S = STATES[cellI];
        const scalarField& A = ALGEBRAIC[cellI];
        const scalarField& R = RATES[cellI];
        // Header line
        Info<< "DEBUG cell=" << cellI
            << " t=" << t1 << "→" << t2;

        if (step >= 0)
            Info<< " step=" << step;

        // Print each selected variable
        forAll(printedNames, k)
        {
            if
            (
                plan.source[k] == ionicModelIO::COL_VM
             || plan.source[k] == ionicModelIO::COL_STATE
            )
            {
                Info<< " " << printedNames[k] << "=" << S[plan.index[k]];
            }
            else if (plan.source[k] == ionicModelIO::COL_RATE)
            {
                Info<< " " << printedNames[k] << "=" << R[plan.index[k]];
            }
            else if (plan.source[k] == ionicModelIO::COL_ALGEBRAIC)
            {
                Info<< " " << printedNames[k] << "=" << A[plan.index[k]];
            }
            else
            {
                FatalErrorInFunction
                    << "Unexpected selection source id " << plan.source[k]
                    << " in debug plan." << exit(FatalError);
            }
        }
        Info<< nl;
    }

    void Foam::ionicModelIO::mapVariableNames
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
        ionicVariableCompatibility::mapVariableNames
        (
            names,
            stateNames,
            nStates,
            algNames,
            nAlg,
            stateIndex,
            algIndex,
            rateIndex
        );
    }







    void Foam::ionicModelIO::writeOneSweepRow
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
        SelectedMapCache& selectedPlanCache
    )
    {
        const SelectionPlan& plan =
            selectedPlan
            (
                deps,
                stateNames,
                nStates,
                algNames,
                nAlg,
                selectedPlanCache
            );

        os << V;

        forAll(plan.source, i)
        {
            if
            (
                plan.source[i] == ionicModelIO::COL_VM
             || plan.source[i] == ionicModelIO::COL_STATE
            )
            {
                os << "," << STATES[plan.index[i]];
            }
            else if (plan.source[i] == ionicModelIO::COL_ALGEBRAIC)
                os << "," << ALG[plan.index[i]];
            else if (plan.source[i] == ionicModelIO::COL_RATE)
            {
                os << "," << RATES[plan.index[i]];
            }
            else
            {
                FatalErrorInFunction
                    << "Unexpected selection source id " << plan.source[i]
                    << " in sweep plan." << exit(FatalError);
            }
        }
        os << nl;
    }
    void Foam::ionicModelIO::writeSweepHeader
    (
        OFstream& os,
        const wordList& deps
    )
    {
        os << "V";
        forAll(deps, i)
        {
            os << "," << deps[i];
        }
        os << nl;
    }



} // End namespace Foam
