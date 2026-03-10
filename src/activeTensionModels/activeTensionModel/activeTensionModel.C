/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  Base class for active tension models
\*---------------------------------------------------------------------------*/

#include "activeTensionModel.H"

namespace Foam
{

defineTypeNameAndDebug(activeTensionModel, 0);
defineRunTimeSelectionTable(activeTensionModel, dictionary);

autoPtr<activeTensionModel> activeTensionModel::New
(
    const dictionary& dict,
    const label nIntegrationPoints
)
{
    const word modelType(dict.lookup("activeTensionModel"));

    Info<< "Selecting activeTensionModel: " << modelType << nl << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown activeTensionModel type: " << modelType << nl
            << "Valid activeTensionModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<activeTensionModel>(cstrIter()(dict, nIntegrationPoints));
}


void activeTensionModel::validateProvider() const
{
    const Requirements req = requirements();

    const bool needsAny = req.needVm || req.needAct || req.needCai;

    if (needsAny && !providerPtr_)
    {
        FatalErrorInFunction
            << "Active tension model requires coupling signals, "
            << "but CouplingSignalProvider was not set."
            << abort(FatalError);
    }

    if (!providerPtr_) return;

    const CouplingSignalProvider& p = *providerPtr_;

    auto require = [&](const CouplingSignal s, const char* name)
    {
        if (!p.hasSignal(s))
        {
            FatalErrorInFunction
                << "Active tension model requires signal '" << name << "' "
                << "but provider does not supply it."
                << abort(FatalError);
        }
    };

    if (req.needVm)  require(CouplingSignal::Vm,  "Vm");
    if (req.needAct) require(CouplingSignal::Act, "Act");
    if (req.needCai) require(CouplingSignal::Cai, "Cai");
}


wordList activeTensionModel::availableFieldNames() const
{
    if (!hasIOMetadata())
        return wordList();

    const label nStates = ioNumStates();
    const label nAlg    = ioNumAlgebraic();
    const char* const* stateNames = ioStateNames();
    const char* const* algNames   = ioAlgebraicNames();

    wordList names(nStates + nAlg + nStates);

    for (label s = 0; s < nStates; ++s)
        names[s] = stateNames[s];

    for (label a = 0; a < nAlg; ++a)
        names[nStates + a] = algNames[a];

    for (label s = 0; s < nStates; ++s)
        names[nStates + nAlg + s] = word("RATES_") + word(stateNames[s]);

    return names;
}


void activeTensionModel::fieldValues
(
    const label integrationPtI,
    scalarField& values
) const
{
    if (!hasIOMetadata())
    {
        values.clear();
        return;
    }

    const auto* statesPtr = ioStatesPtr();
    const auto* algPtr    = ioAlgebraicPtr();
    const auto* ratesPtr  = ioRatesPtr();

    if (!statesPtr || !algPtr || !ratesPtr)
    {
        values.clear();
        return;
    }

    const label nStates = ioNumStates();
    const label nAlg    = ioNumAlgebraic();

    values.setSize(nStates + nAlg + nStates);

    for (label s = 0; s < nStates; ++s)
        values[s] = (*statesPtr)[integrationPtI][s];

    for (label a = 0; a < nAlg; ++a)
        values[nStates + a] = (*algPtr)[integrationPtI][a];

    for (label r = 0; r < nStates; ++r)
        values[nStates + nAlg + r] = (*ratesPtr)[integrationPtI][r];
}


void activeTensionModel::writeHeader(OFstream& os) const
{
    if (!hasIOMetadata())
    {
        FatalErrorInFunction
            << "I/O metadata not available for activeTensionModel "
            << typeName
            << ". Override writeHeader() or provide io* hooks."
            << exit(FatalError);
    }

    const wordList names = exportedFieldNames();

    if (!names.empty())
    {
        activeTensionIO::writeSelectedHeader(os, names);
    }
    else
    {
        activeTensionIO::writeHeader
        (
            os,
            ioStateNames(), ioNumStates(),
            ioAlgebraicNames(), ioNumAlgebraic(),
            fullWritePlanCache_
        );
    }
}


void activeTensionModel::write(const scalar t, OFstream& os) const
{
    if (!hasIOMetadata())
    {
        FatalErrorInFunction
            << "I/O metadata not available for activeTensionModel "
            << typeName
            << ". Override write() or provide io* hooks."
            << exit(FatalError);
    }

    const auto* statesPtr = ioStatesPtr();
    const auto* algPtr    = ioAlgebraicPtr();
    const auto* ratesPtr  = ioRatesPtr();

    if (!statesPtr || !algPtr || !ratesPtr)
    {
        FatalErrorInFunction
            << "I/O storage pointers not available for activeTensionModel "
            << typeName
            << ". Override write() or provide io* hooks."
            << exit(FatalError);
    }

    const wordList names = exportedFieldNames();

    if (!names.empty())
    {
        activeTensionIO::writeSelected
        (
            t, os,
            *statesPtr, *algPtr,
            names,
            ioStateNames(), ioNumStates(),
            ioAlgebraicNames(), ioNumAlgebraic(),
            writeSelectedPlanCache_,
            *ratesPtr
        );
    }
    else
    {
        activeTensionIO::write(t, os, *statesPtr, *algPtr, *ratesPtr);
    }
}


void activeTensionModel::exportStates(PtrList<volScalarField>& outFields)
{
    if (!hasIOMetadata())
    {
        FatalErrorInFunction
            << "I/O metadata not available for activeTensionModel "
            << typeName
            << ". Override exportStates() or provide io* hooks."
            << exit(FatalError);
    }

    const auto* statesPtr = ioStatesPtr();
    const auto* algPtr    = ioAlgebraicPtr();
    const auto* ratesPtr  = ioRatesPtr();

    if (!statesPtr || !algPtr || !ratesPtr)
    {
        FatalErrorInFunction
            << "I/O storage pointers not available for activeTensionModel "
            << typeName
            << ". Override exportStates() or provide io* hooks."
            << exit(FatalError);
    }

    activeTensionIO::exportStateFields
    (
        *statesPtr, *algPtr, *ratesPtr,
        exportedFieldNames(),
        ioStateNames(), ioNumStates(),
        ioAlgebraicNames(), ioNumAlgebraic(),
        exportSelectedPlanCache_,
        outFields
    );
}


void activeTensionModel::debugPrintFields
(
    label cellI,
    scalar t1,
    scalar t2,
    scalar step
) const
{
    if (!hasIOMetadata())
    {
        FatalErrorInFunction
            << "I/O metadata not available for activeTensionModel "
            << typeName
            << ". Override debugPrintFields() or provide io* hooks."
            << exit(FatalError);
    }

    const auto* statesPtr = ioStatesPtr();
    const auto* algPtr    = ioAlgebraicPtr();
    const auto* ratesPtr  = ioRatesPtr();

    if (!statesPtr || !algPtr || !ratesPtr)
    {
        FatalErrorInFunction
            << "I/O storage pointers not available for activeTensionModel "
            << typeName
            << ". Override debugPrintFields() or provide io* hooks."
            << exit(FatalError);
    }

    activeTensionIO::debugPrintFields
    (
        *statesPtr, *algPtr, *ratesPtr,
        debugPrintedNames(),
        ioStateNames(), ioNumStates(),
        ioAlgebraicNames(), ioNumAlgebraic(),
        debugSelectedPlanCache_,
        cellI, t1, t2, step
    );
}

void activeTensionModel::calculateTension
(
    const scalar t,
    const scalar dt,
    const scalarField& lambda,
    scalarField& Ta
)
{
    if (Ta.size() != nIntegrationPoints_)
    {
        FatalErrorInFunction
            << "Ta.size() (" << Ta.size()
            << ") != nIntegrationPoints (" << nIntegrationPoints_ << ")"
            << abort(FatalError);
    }

    currentT_  = t;
    currentDt_ = dt;

    const CouplingSignalProvider& p = provider();
    const CouplingSignal sig = driveSignal();

    forAll(Ta, i)
    {
        currentDriveSignal_ = p.signal(i, sig);
        solveAtPoint(i, currentDriveSignal_, lambda[i], Ta[i]);
    }
}

} // End namespace Foam
