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

#include "electromechanicalVerificationModel.H"

namespace Foam
{

defineTypeNameAndDebug(electromechanicalVerificationModel, 0);
defineRunTimeSelectionTable(electromechanicalVerificationModel, dictionary);


electromechanicalVerificationModel::electromechanicalVerificationModel
(
    const dictionary& dict
)
:
    dict_(dict)
{}


bool electromechanicalVerificationModel::configured(const dictionary& dict)
{
    return dict.found("electromechanicalVerificationModel");
}


autoPtr<electromechanicalVerificationModel>
electromechanicalVerificationModel::New(const dictionary& dict)
{
    const dictionary& verifyDict =
        dict.subDict("electromechanicalVerificationModel");

    const word modelType(verifyDict.lookup("type"));

    Info<< "Selecting electromechanicalVerificationModel: "
        << modelType << nl << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            verifyDict,
            "type",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<electromechanicalVerificationModel>(ctorPtr(verifyDict));
}

} // End namespace Foam

// ************************************************************************* //
