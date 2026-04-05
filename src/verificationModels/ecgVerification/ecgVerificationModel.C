/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "ecgVerificationModel.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(ecgVerificationModel, 0);
defineRunTimeSelectionTable(ecgVerificationModel, dictionary);


ecgVerificationModel::ecgVerificationModel
(
    const electroStateProvider& stateProvider,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
:
    stateProvider_(stateProvider),
    mesh_(stateProvider.mesh()),
    electrodeNames_(electrodeNames),
    electrodePositions_(electrodePositions)
{}


word ecgVerificationModel::selectedType(const dictionary& dict)
{
    const word modelType
    (
        dict.lookupOrDefault<word>("ecgVerificationModel", word::null)
    );

    if (!modelType.empty())
    {
        return modelType;
    }

    if (dict.found("manufactured"))
    {
        return "pseudoECGManufacturedVerifier";
    }

    return word::null;
}


autoPtr<ecgVerificationModel> ecgVerificationModel::New
(
    const electroStateProvider& stateProvider,
    const dictionary& dict,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
{
    const word modelType(selectedType(dict));

    if (modelType.empty())
    {
        return autoPtr<ecgVerificationModel>(nullptr);
    }

    Info<< "Selecting ecgVerificationModel " << modelType << nl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown ecgVerificationModel type " << modelType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    autoPtr<ecgVerificationModel> verifier
    (
        ctorPtr(stateProvider, dict, electrodeNames, electrodePositions)
    );
    verifier->validateProvider();
    return verifier;
}


const volScalarField& ecgVerificationModel::requireVm() const
{
    const volScalarField* fieldPtr = stateProvider_.VmPtr();

    if (!fieldPtr)
    {
        FatalErrorInFunction
            << "ECG verification model requires Vm, but the selected "
            << "electroStateProvider does not expose it."
            << exit(FatalError);
    }

    return *fieldPtr;
}


const volScalarField& ecgVerificationModel::requirePhiE() const
{
    const volScalarField* fieldPtr = stateProvider_.phiEPtr();

    if (!fieldPtr)
    {
        FatalErrorInFunction
            << "ECG verification model requires phiE, but the selected "
            << "electroStateProvider does not expose it."
            << exit(FatalError);
    }

    return *fieldPtr;
}


const volTensorField& ecgVerificationModel::requireConductivity() const
{
    const volTensorField* fieldPtr = stateProvider_.conductivityPtr();

    if (!fieldPtr)
    {
        FatalErrorInFunction
            << "ECG verification model requires conductivity, but the selected "
            << "electroStateProvider does not expose it."
            << exit(FatalError);
    }

    return *fieldPtr;
}


void ecgVerificationModel::validateProvider() const
{
    const Requirements needs = requirements();

    if (needs.needVm)
    {
        (void)requireVm();
    }

    if (needs.needPhiE)
    {
        (void)requirePhiE();
    }

    if (needs.needConductivity)
    {
        (void)requireConductivity();
    }
}


void ecgVerificationModel::updateElectrodes
(
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
{
    electrodeNames_ = electrodeNames;
    electrodePositions_ = electrodePositions;
}

} // End namespace Foam

// ************************************************************************* //
