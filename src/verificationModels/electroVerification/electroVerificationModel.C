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

#include "electroVerificationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

namespace Foam
{

defineTypeNameAndDebug(electroVerificationModel, 0);
defineRunTimeSelectionTable(electroVerificationModel, dictionary);

electroVerificationModel::electroVerificationModel(const dictionary& dict)
:
    dict_(dict)
{}

autoPtr<electroVerificationModel> electroVerificationModel::New
(
    const dictionary& dict
)
{
    const word modelType(selectedType(dict));

    if (modelType.empty() || modelType == "none")
    {
        return autoPtr<electroVerificationModel>(nullptr);
    }

    Info<< "Selecting electroVerificationModel " << modelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown electroVerificationModel type " << modelType
            << nl << nl << "Valid electroVerificationModels are : " << nl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<electroVerificationModel>(cstrIter()(dict));
}

word electroVerificationModel::selectedType(const dictionary& dict)
{
    const dictionary* subDictPtr = dict.findDict("verificationModel");
    if (subDictPtr)
    {
        return subDictPtr->lookupOrDefault<word>("type", "");
    }
    return "";
}

void electroVerificationModel::allocateFields
(
    const wordList& fieldNames,
    const fvMesh& mesh,
    const word& prefix,
    PtrList<volScalarField>& fields
)
{
    fields.setSize(fieldNames.size());
    forAll(fieldNames, i)
    {
        fields.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    prefix + fieldNames[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimless,
                "zeroGradient"
            )
        );
    }
}

} // End namespace Foam

// ************************************************************************* //
