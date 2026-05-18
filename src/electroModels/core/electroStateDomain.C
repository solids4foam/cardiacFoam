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

#include "electroStateDomain.H"
#include "error.H"

namespace Foam
{

defineTypeNameAndDebug(electroStateDomain, 0);
defineRunTimeSelectionTable(electroStateDomain, dictionary);


autoPtr<electroStateDomain> electroStateDomain::New
(
    const fvMesh& baseMesh,
    myocardiumDomainInterface& myocardium,
    const dictionary& dict
)
{
    const word domainType(dict.lookup("type"));

    Info<< nl << "Selecting electroStateDomain " << domainType << endl;

    auto* ctorPtr = dictionaryConstructorTable(domainType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "type",
            domainType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<electroStateDomain>(ctorPtr(baseMesh, myocardium, dict));
}

} // End namespace Foam

// ************************************************************************* //
