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

#include "eikonalVerificationModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(eikonalVerificationModel, 0);
defineRunTimeSelectionTable(eikonalVerificationModel, dictionary);

eikonalVerificationModel::eikonalVerificationModel(const dictionary& dict)
:
    dict_(dict)
{}


word eikonalVerificationModel::selectedType(const dictionary& dict)
{
    const dictionary* subDictPtr = dict.findDict("verificationModel");
    if (subDictPtr)
    {
        return subDictPtr->lookupOrDefault<word>("type", "");
    }
    return "";
}


autoPtr<eikonalVerificationModel> eikonalVerificationModel::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volTensorField& conductivity,
    const dimensionedScalar& chi,
    const dimensionedScalar& Cm,
    const dimensionedScalar& c0,
    const Switch& eikonalAdvectionDiffusionApproach
)
{
    const word modelType(selectedType(dict));

    if (modelType.empty() || modelType == "none")
    {
        return autoPtr<eikonalVerificationModel>(nullptr);
    }

    Info<< "Selecting eikonalVerificationModel " << modelType << endl;

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "eikonalVerificationModel table is empty. "
            << "Did you forget to load libverificationModels.so in controlDict?"
            << exit(FatalError);
    }

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown eikonalVerificationModel type " << modelType
            << nl << nl << "Valid eikonalVerificationModels are : " << nl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<eikonalVerificationModel>
    (
        cstrIter()
        (
            dict,
            mesh,
            conductivity,
            chi,
            Cm,
            c0,
            eikonalAdvectionDiffusionApproach
        )
    );
}

} // End namespace Foam

// ************************************************************************* //
