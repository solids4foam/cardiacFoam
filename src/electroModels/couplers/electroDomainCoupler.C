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

#include "electroDomainCoupler.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(ElectroDomainCoupler, 0);
defineRunTimeSelectionTable(ElectroDomainCoupler, dictionary);


autoPtr<ElectroDomainCoupler> ElectroDomainCoupler::New
(
    tissueCouplingEndpoint& primaryDomain,
    electrophysicsDomain&       secondaryDomain,
    const dictionary&         dict
)
{
    const word modelType
    (
        dict.lookupOrDefault<word>("ElectroDomainCoupler", word::null)
    );

    if (modelType.empty())
    {
        FatalErrorInFunction
            << "No electro-domain coupling model specified in dictionary '"
            << dict.dictName() << "'." << nl
            << "Expected key ElectroDomainCoupler."
            << exit(FatalError);
    }

    Info<< "Selecting ElectroDomainCoupler " << modelType << nl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown ElectroDomainCoupler type " << modelType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ElectroDomainCoupler>
    (
        ctorPtr(primaryDomain, secondaryDomain, dict)
    );
}

} // End namespace Foam

// ************************************************************************* //
