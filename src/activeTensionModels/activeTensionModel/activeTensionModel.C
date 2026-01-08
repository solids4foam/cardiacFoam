
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

} // namespace Foam
