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

#include "conductionSystemDomain.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(ConductionSystemDomain, 0);
defineRunTimeSelectionTable(ConductionSystemDomain, dictionary);


ConductionSystemDomain::ConductionSystemDomain
(
    const ConductionSystemDomainContext& context,
    const dictionary&     dict
)
:
    context_(context)
{}


autoPtr<ConductionSystemDomain>
ConductionSystemDomain::New
(
    const ConductionSystemDomainContext& context,
    const dictionary&                    dict,
    scalar                               initialDeltaT
)
{
    const word modelType
    (
        dict.lookupOrDefault<word>("ConductionSystemDomain", word::null)
    );

    if (modelType.empty())
    {
        FatalErrorInFunction
            << "No ConductionSystemDomain specified in dictionary '"
            << dict.dictName() << "'." << nl
            << "Expected key ConductionSystemDomain."
            << exit(FatalError);
    }

    Info<< "Selecting ConductionSystemDomain " << modelType
        << nl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown ConductionSystemDomain type "
            << modelType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ConductionSystemDomain>
    (
        ctorPtr(context, dict, initialDeltaT)
    );
}

} // End namespace Foam

// ************************************************************************* //
