/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  I/O helpers for active tension models (mirrors ionicModelIO pattern).
\*---------------------------------------------------------------------------*/

#include "activeTensionIO.H"

namespace Foam
{

namespace
{

// Emit header line: "time col1 col2 ..."
void emitHeader(OFstream& os, const wordList& names)
{
    os << "time";
    forAll(names, k)
    {
        os << " " << names[k];
    }
    os << nl;
}

// Emit one data row using a resolved selection plan
void emitRow
(
    const scalar t,
    OFstream& os,
    const scalarField& S,
    const scalarField& A,
    const scalarField& R,
    const activeTensionIO::SelectionPlan& plan
)
{
    os << t;
    forAll(plan.source, k)
    {
        const label src = plan.source[k];
        const label idx = plan.index[k];

        if (src == activeTensionIO::COL_STATE)
        {
            os << " " << S[idx];
        }
        else if (src == activeTensionIO::COL_ALGEBRAIC)
        {
            os << " " << A[idx];
        }
        else if (src == activeTensionIO::COL_RATE)
        {
            os << " " << R[idx];
        }
        else
        {
            FatalErrorInFunction
                << "Unexpected column source id " << src
                << " in active tension write plan."
                << exit(FatalError);
        }
    }
    os << nl;
}

// Build or retrieve a cached selection plan for a given set of exported names
const activeTensionIO::SelectionPlan& buildSelectedPlan
(
    const wordList& exportedNames,
    const char* const stateNames[],
    const label nStates,
    const char* const algNames[],
    const label nAlg,
    activeTensionIO::SelectedMapCache& cache
)
{
    const bool cacheHit =
           cache.valid
        && cache.stateNamesPtr == static_cast<const void*>(stateNames)
        && cache.algNamesPtr   == static_cast<const void*>(algNames)
        && cache.nStates == nStates
        && cache.nAlg    == nAlg
        && cache.exportedNames == exportedNames;

    if (!cacheHit)
    {
        cache.stateNamesPtr = static_cast<const void*>(stateNames);
        cache.algNamesPtr   = static_cast<const void*>(algNames);
        cache.nStates = nStates;
        cache.nAlg    = nAlg;
        cache.exportedNames = exportedNames;

        activeTensionIO::SelectionPlan& plan = cache.plan;
        plan.names = exportedNames;
        plan.source.setSize(exportedNames.size());
        plan.index.setSize(exportedNames.size());

        const word ratePrefix("RATES_");

        forAll(exportedNames, k)
        {
            const word& name = exportedNames[k];
            label stateIdx = -1;
            label algIdx   = -1;
            label rateIdx  = -1;

            for (label i = 0; i < nStates && stateIdx < 0; ++i)
            {
                if (word(stateNames[i]) == name)
                    stateIdx = i;
            }
            for (label i = 0; i < nAlg && stateIdx < 0 && algIdx < 0; ++i)
            {
                if (word(algNames[i]) == name)
                    algIdx = i;
            }
            if (stateIdx < 0 && algIdx < 0)
            {
                for (label i = 0; i < nStates; ++i)
                {
                    if ((ratePrefix + word(stateNames[i])) == name)
                    {
                        rateIdx = i;
                        break;
                    }
                }
            }

            if (stateIdx >= 0)
            {
                plan.source[k] = activeTensionIO::COL_STATE;
                plan.index[k]  = stateIdx;
            }
            else if (algIdx >= 0)
            {
                plan.source[k] = activeTensionIO::COL_ALGEBRAIC;
                plan.index[k]  = algIdx;
            }
            else if (rateIdx >= 0)
            {
                plan.source[k] = activeTensionIO::COL_RATE;
                plan.index[k]  = rateIdx;
            }
            else
            {
                WarningInFunction
                    << "Unknown active tension variable '" << name
                    << "'. Output will be zero." << nl;
                plan.source[k] = activeTensionIO::COL_STATE;
                plan.index[k]  = 0;
            }
        }

        cache.valid = true;
    }

    return cache.plan;
}

// Build or retrieve a cached full plan (all states + algebraic + rates)
const activeTensionIO::SelectionPlan& buildFullPlan
(
    const char* const stateNames[],
    const label nStates,
    const char* const algNames[],
    const label nAlg,
    activeTensionIO::FullPlanCache& cache
)
{
    const bool cacheHit =
           cache.valid
        && cache.stateNamesPtr == static_cast<const void*>(stateNames)
        && cache.algNamesPtr   == static_cast<const void*>(algNames)
        && cache.nStates == nStates
        && cache.nAlg    == nAlg;

    if (!cacheHit)
    {
        cache.stateNamesPtr = static_cast<const void*>(stateNames);
        cache.algNamesPtr   = static_cast<const void*>(algNames);
        cache.nStates = nStates;
        cache.nAlg    = nAlg;

        activeTensionIO::SelectionPlan& plan = cache.plan;
        const label nCols = nStates + nAlg + nStates;
        plan.names.setSize(nCols);
        plan.source.setSize(nCols);
        plan.index.setSize(nCols);

        label k = 0;

        for (label s = 0; s < nStates; ++s)
        {
            plan.names[k]  = stateNames[s];
            plan.source[k] = activeTensionIO::COL_STATE;
            plan.index[k]  = s;
            ++k;
        }

        for (label a = 0; a < nAlg; ++a)
        {
            plan.names[k]  = algNames[a];
            plan.source[k] = activeTensionIO::COL_ALGEBRAIC;
            plan.index[k]  = a;
            ++k;
        }

        for (label s = 0; s < nStates; ++s)
        {
            plan.names[k]  = word("RATES_") + word(stateNames[s]);
            plan.source[k] = activeTensionIO::COL_RATE;
            plan.index[k]  = s;
            ++k;
        }

        cache.valid = true;
    }

    return cache.plan;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

void activeTensionIO::writeHeader
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
        buildFullPlan(stateNames, nStates, algNames, nAlg, fullPlanCache).names
    );
}


void activeTensionIO::writeSelectedHeader
(
    OFstream& os,
    const wordList& exportedNames
)
{
    emitHeader(os, exportedNames);
}


void activeTensionIO::write
(
    scalar t,
    OFstream& os,
    const PtrList<scalarField>& STATES,
    const PtrList<scalarField>& ALGEBRAIC,
    const PtrList<scalarField>& RATES
)
{
    const scalarField& S = STATES[0];
    const scalarField& A = ALGEBRAIC[0];
    const scalarField& R = RATES[0];

    os << t;
    forAll(S, i) { os << " " << S[i]; }
    forAll(A, i) { os << " " << A[i]; }
    forAll(R, i) { os << " " << R[i]; }
    os << nl;
}


void activeTensionIO::writeSelected
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
    const PtrList<scalarField>& RATES
)
{
    const SelectionPlan& plan =
        buildSelectedPlan
        (
            exportedNames,
            stateNames, nStates,
            algNames,   nAlg,
            selectedPlanCache
        );

    emitRow(t, os, STATES[0], ALGEBRAIC[0], RATES[0], plan);
}


const Foam::wordList& activeTensionIO::exportedFieldNamesRef
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
        && cache.algNamesPtr   == static_cast<const void*>(algNames)
        && cache.nStates == nStates
        && cache.nAlg    == nAlg
        && cache.requested == userList;

    if (cacheHit)
        return cache.filtered;

    const word ratePrefix("RATES_");
    wordList filtered;

    forAll(userList, k)
    {
        const word& name = userList[k];
        bool found = false;

        for (label i = 0; i < nStates && !found; ++i)
            if (word(stateNames[i]) == name) found = true;

        for (label i = 0; i < nAlg && !found; ++i)
            if (word(algNames[i]) == name) found = true;

        for (label i = 0; i < nStates && !found; ++i)
            if ((ratePrefix + word(stateNames[i])) == name) found = true;

        if (found)
            filtered.append(name);
        else
            WarningInFunction
                << "Ignoring unknown active tension variable '" << name
                << "'." << nl;
    }

    cache.valid         = true;
    cache.stateNamesPtr = static_cast<const void*>(stateNames);
    cache.algNamesPtr   = static_cast<const void*>(algNames);
    cache.nStates       = nStates;
    cache.nAlg          = nAlg;
    cache.requested     = userList;
    cache.filtered      = filtered;

    return cache.filtered;
}


void activeTensionIO::exportStateFields
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
        buildSelectedPlan
        (
            exportedNames,
            stateNames, nStates,
            algNames,   nAlg,
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

    forAll(STATES, cellI)
    {
        const scalarField& S = STATES[cellI];
        const scalarField& A = ALGEBRAIC[cellI];
        const scalarField& R = RATES[cellI];

        forAll(outFields, k)
        {
            const label src = plan.source[k];
            const label idx = plan.index[k];

            if (src == COL_STATE)
                outFields[k][cellI] = S[idx];
            else if (src == COL_ALGEBRAIC)
                outFields[k][cellI] = A[idx];
            else if (src == COL_RATE)
                outFields[k][cellI] = R[idx];
            else
            {
                FatalErrorInFunction
                    << "Unexpected column source id " << src
                    << " in active tension export plan."
                    << exit(FatalError);
            }
        }
    }

    forAll(outFields, k)
    {
        outFields[k].correctBoundaryConditions();
    }
}


void activeTensionIO::debugPrintFields
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
        buildSelectedPlan
        (
            printedNames,
            stateNames, nStates,
            algNames,   nAlg,
            selectedPlanCache
        );

    const scalarField& S = STATES[cellI];
    const scalarField& A = ALGEBRAIC[cellI];
    const scalarField& R = RATES[cellI];

    Info<< "DEBUG cell=" << cellI
        << " t=" << t1 << "->" << t2;

    if (step >= 0)
        Info<< " step=" << step;

    forAll(printedNames, k)
    {
        const label src = plan.source[k];
        const label idx = plan.index[k];

        if (src == COL_STATE)
            Info<< " " << printedNames[k] << "=" << S[idx];
        else if (src == COL_ALGEBRAIC)
            Info<< " " << printedNames[k] << "=" << A[idx];
        else if (src == COL_RATE)
            Info<< " " << printedNames[k] << "=" << R[idx];
    }

    Info<< nl;
}

} // End namespace Foam
