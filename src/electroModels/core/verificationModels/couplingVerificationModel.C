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
    const word modelType(selectedType(dict));

    Info<< "Selecting couplingVerificationModel " << modelType << nl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);
    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown couplingVerificationModel type " << modelType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<couplingVerificationModel>(cstrIter()(dict));
}


word couplingVerificationModel::selectedType(const dictionary& dict)
{
    return dict.get<word>("type");
}


couplingVerificationModel::couplingVerificationModel(const dictionary& dict)
:
    dict_(dict)
{}

} // End namespace Foam

// ************************************************************************* //
