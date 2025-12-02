#include "ionicModelIO.H"

namespace Foam {

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
        if (outFields.size() != exportedNames.size())
        {
            FatalErrorInFunction
                << "exportStateFields: expected " << exportedNames.size()
                << " fields, got " << outFields.size()
                << exit(FatalError);
        }

        // Precompute mapping name → state or algebraic index
        List<label> stateIndex(exportedNames.size(), -1);
        List<label> algIndex(exportedNames.size(), -1);

        forAll(exportedNames, k)
        {
            const word& name = exportedNames[k];

            // Match STATE variables
            for (label s = 0; s < nStates; ++s)
            {
                if (name == stateNames[s])
                {
                    stateIndex[k] = s;
                    goto mapped;
                }
            }

            // Match ALGEBRAIC variables
            for (label a = 0; a < nAlg; ++a)
            {
                if (name == algNames[a])
                {
                    algIndex[k] = a;
                    goto mapped;
                }
            }

            FatalErrorInFunction
                << "Unknown variable '" << name
                << "' in exportedVariables list"
                << exit(FatalError);

        mapped:;
        }

        // Populate volScalarFields
        forAll(STATES, cellI)
        {
            const scalarField& S = STATES[cellI];
            const scalarField& A = ALGEBRAIC[cellI];

            forAll(outFields, k)
            {
                if (stateIndex[k] >= 0)
                {
                    (*outFields[k])[cellI] = S[stateIndex[k]];
                }
                else
                {
                    (*outFields[k])[cellI] = A[algIndex[k]];
                }
            }
        }

        // Fix boundaries
        forAll(outFields, k)
            outFields[k]->correctBoundaryConditions();
    }

}



    