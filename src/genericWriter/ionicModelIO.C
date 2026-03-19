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

namespace Foam {

    namespace
    {
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
        emitHeader(os, fullPlan(stateNames, nStates, algNames, nAlg, fullPlanCache).names);
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
        return (tEnd >= writeAfterTime);
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

        forAll(STATES, cellI)
        {
            scalarField& S = STATES[cellI];

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
            if (plan.source[k] == ionicModelIO::COL_VM || plan.source[k] == ionicModelIO::COL_STATE)
                {Info<< " " << printedNames[k] << "=" << S[plan.index[k]];}
            else if (plan.source[k] == ionicModelIO::COL_RATE)
            {
                Info<< " " << printedNames[k] << "=" << R[plan.index[k]];
            }
            else if (plan.source[k] == ionicModelIO::COL_ALGEBRAIC)
                {Info<< " " << printedNames[k] << "=" << A[plan.index[k]];}
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
            if (plan.source[i] == ionicModelIO::COL_VM || plan.source[i] == ionicModelIO::COL_STATE)
                os << "," << STATES[plan.index[i]];
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







