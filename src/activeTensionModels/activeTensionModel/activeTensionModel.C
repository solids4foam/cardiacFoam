
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
#include "activeTensionIO.H"

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

    const bool needsAny =
        req.needVm || req.needAct || req.needCai;

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

    if (req.needVm)       require(CouplingSignal::Vm,       "Vm");
    if (req.needAct)      require(CouplingSignal::Act,      "Act");
    if (req.needCai)      require(CouplingSignal::Cai,      "Cai");
}

void activeTensionModel::writeHeader(OFstream& os) const
{
    wordList names = exportedFieldNames();
    if (names.empty())
    {
        names = availableFieldNames();
    }

    if (!names.empty())
    {
        activeTensionIO::writeHeader(os, names);
    }
    else
    {
        os << "# t Ta" << nl;
    }
}

void activeTensionModel::write(const scalar t, OFstream& os) const
{
    scalarField values;
    fieldValues(0, values);

    if (values.empty())
    {
        os << t << " " << 0.0 << nl;
        return;
    }

    wordList names = exportedFieldNames();
    if (names.empty())
    {
        activeTensionIO::write(t, os, values);
        return;
    }

    wordList available = availableFieldNames();
    scalarField selected(names.size(), 0.0);
    forAll(names, i)
    {
        label idx = available.find(names[i]);
        if (idx >= 0 && idx < values.size())
        {
            selected[i] = values[idx];
        }
    }
    activeTensionIO::write(t, os, selected);
}

} // namespace Foam
