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

#include "graphVerificationModel.H"

namespace Foam
{

defineRunTimeSelectionTable(graphVerificationModel, dictionary);
defineTypeNameAndDebug(graphVerificationModel, 0);

autoPtr<graphVerificationModel> graphVerificationModel::New(const dictionary& dict)
{
    const word modelType(selectedType(dict));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown graphVerificationModel type "
            << modelType << nl << nl
            << "Valid graphVerificationModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<graphVerificationModel>(cstrIter()(dict));
}


word graphVerificationModel::selectedType(const dictionary& dict)
{
    return dict.get<word>("type");
}


graphVerificationModel::graphVerificationModel(const dictionary& dict)
:
    dict_(dict)
{}

} // End namespace Foam

// ************************************************************************* //
