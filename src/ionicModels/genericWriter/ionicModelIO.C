#include "ionicModelIO.H"

namespace Foam {

    fileName ionicModelIO::createOutputFile
    (
        const dictionary& dict,
        const word& modelName,
        const word& tissueName,
        const Time& runTime
    )
    {
        label s1 = readLabel(dict.lookup("stim_period_S1"));
        label s2 = dict.found("stim_period_S2")
                    ? readLabel(dict.lookup("stim_period_S2"))
                    : -1;

        const bool protocolMode = (s2 > 0);

        if (protocolMode)
        {
            return runTime.path()
                / (modelName + "_" + tissueName + "_"
                + name(s1) + "-" + name(s2) + "ms.txt");
        }
        else
        {
            return runTime.path()
                / (modelName + "_" + tissueName + ".txt");
        }
    }

    bool ionicModelIO::shouldWriteStep
    (
        const dictionary& dict,
        const scalar currentTime
    )
    {
        label s2 = dict.found("stim_period_S2")
            ? readLabel(dict.lookup("stim_period_S2"))
            : -1;

        const bool protocolMode = (s2 > 0);

        if (!protocolMode)
        {
            // Regular pacing → always write
            return true;
        }

        // S1–S2 protocol: only write after 8 seconds
        return currentTime > 8.0;
    }

    void ionicModelIO::writeTimestep
    (
        ionicModel& model,
        OFstream& os,
        const dictionary& dict,
        const Time& runTime
    )
    {
        if (!shouldWriteStep(dict, runTime.value()))
        {
            return;
        }

        model.write(runTime.value(), os);
    }


    void ionicModelIO::loadStimulusConstants
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
    )
    {
        const char* requiredKeys[] =
        {
            "stim_start",
            "stim_period_S1",
            "stim_duration",
            "stim_amplitude"
        };

        for (const char* k : requiredKeys)
        {
            if (!dict.found(k))
            {
                FatalErrorInFunction
                    << "Missing required stimulus key: " << k
                    << exit(FatalError);
            }
        }

        CONSTANTS[stim_start]      = readScalar(dict.lookup("stim_start"));
        CONSTANTS[stim_period_S1]  = readScalar(dict.lookup("stim_period_S1"));
        CONSTANTS[stim_duration]   = readScalar(dict.lookup("stim_duration"));
        CONSTANTS[stim_amplitude]  = readScalar(dict.lookup("stim_amplitude"));

        // Optional keys
        if (dict.found("nstim1"))
            CONSTANTS[nstim1] = readScalar(dict.lookup("nstim1"));
        if (dict.found("stim_period_S2"))
            CONSTANTS[stim_period_S2] = readScalar(dict.lookup("stim_period_S2"));
        if (dict.found("nstim2"))
            CONSTANTS[nstim2] = readScalar(dict.lookup("nstim2"));
    }



    void ionicModelIO::writeHeader
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


    void ionicModelIO::write
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

    Foam::wordList Foam::ionicModelIO::exportedFieldNames
    (
    const wordList& userList,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg
    )
    {
        if (userList.size())
        {
            return userList;
        }
        return Foam::wordList();
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
        List<volScalarField*>& outFields
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
                    (*outFields[k])[cellI] = S[stateIndex[k]];
                else
                    (*outFields[k])[cellI] = A[algIndex[k]];
            }
        }

        // 3. Boundaries
        forAll(outFields, k)
            outFields[k]->correctBoundaryConditions();
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
}





