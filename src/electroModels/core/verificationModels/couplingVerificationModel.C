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

#include "couplingVerificationModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(couplingVerificationModel, 0);
defineRunTimeSelectionTable(couplingVerificationModel, dictionary);


autoPtr<couplingVerificationModel> couplingVerificationModel::New
(
    const dictionary& dict
)
{
    const word modelType(dict.lookupOrDefault<word>("type", word::null));

    if (modelType.empty())
    {
        FatalErrorInFunction
            << "No verification model specified in dictionary '"
            << dict.dictName() << "'." << nl
            << "Expected key 'type'."
            << exit(FatalError);
    }

    Info<< "Selecting couplingVerificationModel " << modelType << nl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown couplingVerificationModel type " << modelType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<couplingVerificationModel>(ctorPtr(dict));
}


word couplingVerificationModel::selectedType(const dictionary& dict)
{
    return dict.lookupOrDefault<word>("type", word::null);
}


couplingVerificationModel::couplingVerificationModel(const dictionary& dict)
:
    dict_(dict)
{}

} // End namespace Foam

// ************************************************************************* //
