#include "ionicModelIO.H"

namespace Foam {

    void Foam::ionicModelIO::writeHeader
    (
        OFstream& os,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg
    )
    {
        os << "time Vm";
        // States
        for (int i = 1; i < nStates; ++i)
        {
            os << " " << stateNames[i];
        }
        // Algebraic
        for (int i = 0; i < nAlg; ++i)
        {
            os << " " << algNames[i];
        }
        // Rates
        for (int i = 0; i < nStates; ++i)
        {
            os << " RATES_" << stateNames[i];
        }

        os << endl;
    }

    void Foam::ionicModelIO::writeSelectedHeader
    (
        OFstream& os,
        const wordList& exportedNames
    )
    {
        os << "time";

        forAll(exportedNames, k)
        {
            os << " " << exportedNames[k];
        }

        os << nl;
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

        // time + transformed Vm
        os << t << " " << transformVm(S);

        // states
        for (int i = 1; i < S.size(); ++i)
        {
            os << " " << S[i];
        }

        // algebraic
        const scalarField& A = ALGEBRAIC[0];
        forAll(A, i)
        {
            os << " " << A[i];
        }

        // rates
        const scalarField& R = RATES[0];
        forAll(R, i)
        {
            os << " " << R[i];
        }

        os << endl;
    }

    void Foam::ionicModelIO::writeSelected
    (
        const scalar t,
        OFstream& os,
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const wordList& exportedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg
    )
    {
        // 1. mapping (same as exportStateFields)
        List<label> stateIndex, algIndex;

        mapVariableNames
        (
            exportedNames,
            stateNames, nStates,
            algNames,  nAlg,
            stateIndex,
            algIndex
        );

        // 2. single cell only
        const scalarField& S = STATES[0];
        const scalarField& A = ALGEBRAIC[0];

        // 3. write
        os << t;

        forAll(stateIndex, k)
        {
            if (stateIndex[k] >= 0)
            {
                os << " " << S[stateIndex[k]];
            }
            else
            {
                os << " " << A[algIndex[k]];
            }
        }

        os << nl;
    }


    Foam::wordList Foam::ionicModelIO::exportedFieldNames
    (
        const wordList& userList,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg
    )
    {
        if (!userList.size())
        {
            return Foam::wordList();
        }

        wordList filtered;
        forAll(userList, k)
        {
            const word& name = userList[k];
            bool found = false;
            for (label i = 0; i < nStates; ++i)
            {
                if (name == stateNames[i])
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                for (label i = 0; i < nAlg; ++i)
                {
                    if (name == algNames[i])
                    {
                        found = true;
                        break;
                    }
                }
            }
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

        return filtered;
    }


    void Foam::ionicModelIO::exportStateFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const wordList& exportedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        PtrList<volScalarField>& outFields
    )
    {
        // 1. mapping
        List<label> stateIndex, algIndex;
        mapVariableNames
        (
            exportedNames,
            stateNames, nStates,
            algNames,  nAlg,
            stateIndex,
            algIndex
        );

        // 2. Populate volScalarFields
        forAll(STATES, cellI)
        {
            const scalarField& S = STATES[cellI];
            const scalarField& A = ALGEBRAIC[cellI];

            forAll(outFields, k)
            {
                if (stateIndex[k] >= 0)
                {
                    outFields[k][cellI] = S[stateIndex[k]];
                }
                else
                {
                    outFields[k][cellI] = A[algIndex[k]];
                }
            }
        }

        // 3. Boundaries
        forAll(outFields, k)
        {
            outFields[k].correctBoundaryConditions();
        }
    }

    void Foam::ionicModelIO::debugPrintFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const wordList& printedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        label cellI,
        scalar t1,
        scalar t2,
        scalar step
    )
    {
        if (printedNames.empty())
            return;

        // Map names → indices
        List<label> stateIndex, algIndex;
        mapVariableNames
        (
            printedNames,
            stateNames, nStates,
            algNames,  nAlg,
            stateIndex,
            algIndex
        );

        const scalarField& S = STATES[cellI];
        const scalarField& A = ALGEBRAIC[cellI];
        // Header line
        Info<< "DEBUG cell=" << cellI
            << " t=" << t1 << "→" << t2;

        if (step >= 0)
            Info<< " step=" << step;

        // Print each selected variable
        forAll(printedNames, k)
        {
            if (stateIndex[k] >= 0)
                {Info<< " " << printedNames[k] << "=" << S[stateIndex[k]];}
            else
                {Info<< " " << printedNames[k] << "=" << A[algIndex[k]];}
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
        List<label>& algIndex
    )
    {
        stateIndex.setSize(names.size(), -1);
        algIndex.setSize(names.size(), -1);

        forAll(names, k)
        {
            const word& name = names[k];

            // Match STATE variable
            for (label s = 0; s < nStates; ++s)
            {
                if (name == stateNames[s])
                {
                    stateIndex[k] = s;
                    goto mapped;
                }
            }

            // Match ALGEBRAIC variable
            for (label a = 0; a < nAlg; ++a)
            {
                if (name == algNames[a])
                {
                    algIndex[k] = a;
                    goto mapped;
                }
            }

            FatalErrorInFunction
                << "Unknown debug variable: " << name
                << exit(FatalError);

        mapped:;
        }
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
        label nAlg
    )
    {
        List<label> stateIndex, algIndex;
    
        mapVariableNames(
            deps,
            stateNames, nStates,
            algNames, nAlg,
            stateIndex, algIndex
        );
    
        os << V;
    
        forAll(deps, i)
        {
            if (stateIndex[i] >= 0)
                os << "," << STATES[stateIndex[i]];
            else if (algIndex[i] >= 0)
                os << "," << ALG[algIndex[i]];
            else
                os << ",0";
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






    
