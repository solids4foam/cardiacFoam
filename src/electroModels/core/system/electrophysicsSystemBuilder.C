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

#include "electrophysicsSystemBuilder.H"
#include "conductionSystemDomain.H"
#include "ecgDomain.H"
#include "electroDomainCoupler.H"
#include "electrophysicsAdvanceScheme.H"
#include "error.H"

#include "DynamicList.H"

namespace Foam
{
namespace electrophysicsSystemBuilder
{

// Terminology used by the orchestration layer:
//   primary domain        = myocardium
//   conduction domains    = pre-myocardium graph / Purkinje domains
//   ECG domains           = post-myocardium ECG domains

namespace
{
bool appendSubDictionaries
(
    const dictionary& parent,
    const word& entryName,
    DynamicList<word>& names,
    DynamicList<const dictionary*>& dicts
)
{
    if (!parent.found(entryName))
    {
        return false;
    }

    const dictionary& container = parent.subDict(entryName);
    const label initialSize = names.size();

    forAllConstIter(dictionary, container, iter)
    {
        const entry& e = iter();

        if (!e.isDict())
        {
            continue;
        }

        names.append(e.keyword());
        dicts.append(&e.dict());
    }

    return names.size() > initialSize;
}


void collectConductionDomainDicts
(
    const dictionary& electroProperties,
    DynamicList<word>& names,
    DynamicList<const dictionary*>& dicts
)
{
    appendSubDictionaries
    (
        electroProperties,
        "conductionNetworkDomains",
        names,
        dicts
    );
}


void collectConductionCouplingDicts
(
    const dictionary& electroProperties,
    DynamicList<const dictionary*>& dicts
)
{
    DynamicList<word> ignoredNames;
    appendSubDictionaries
    (
        electroProperties,
        "domainCouplings",
        ignoredNames,
        dicts
    );
}

} // End anonymous namespace


void configureMyocardiumDomain
(
    electrophysicsSystem&    system,
    const fvMesh&            mesh,
    const dictionary&        electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList&          postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    autoPtr<ionicModel>&     ionicModelPtr,
    autoPtr<electroVerificationModel>& verificationModelPtr,
    scalar                   initialDeltaT
)
{
    system.setMyocardium
    (
        myocardiumDomainInterface::New
        (
            mesh,
            electroProperties,
            outFields,
            postProcessFieldNames,
            postProcessFields,
            ionicModelPtr,
            verificationModelPtr,
            initialDeltaT
        ).ptr()
    );
}


void configureAdvanceScheme
(
    electrophysicsSystem& system,
    const dictionary&     electroProperties
)
{
    system.setAdvanceScheme
    (
        electrophysicsAdvanceScheme::New(electroProperties).ptr()
    );
}


void configureConductionDomains
(
    electrophysicsSystem& system,
    const fvMesh&         mesh,
    const dictionary&     electroProperties,
    scalar                initialDeltaT
)
{
    system.clearConductionDomains();
    system.clearConductionCouplings();

    DynamicList<word> conductionDomainNames;
    DynamicList<const dictionary*> conductionDomainDicts;
    collectConductionDomainDicts
    (
        electroProperties,
        conductionDomainNames,
        conductionDomainDicts
    );

    if (conductionDomainNames.empty())
    {
        return;
    }

    HashTable<ConductionSystemDomain*> conductionDomainsByName
    (
        conductionDomainNames.size()
    );

    forAll(conductionDomainNames, i)
    {
        const word& domainName = conductionDomainNames[i];

        if (conductionDomainsByName.found(domainName))
        {
            FatalErrorInFunction
                << "Duplicate conduction domain name '" << domainName
                << "' while configuring pre-myocardium domains."
                << exit(FatalError);
        }

        autoPtr<ConductionSystemDomain> conductionDomain
        (
            ConductionSystemDomain::New
            (
                mesh,
                *conductionDomainDicts[i],
                initialDeltaT
            )
        );

        ConductionSystemDomain* domainPtr = conductionDomain.ptr();
        conductionDomainsByName.insert(domainName, domainPtr);
        system.appendConductionDomain(domainPtr);
    }

    DynamicList<const dictionary*> couplingDicts;
    collectConductionCouplingDicts
    (
        electroProperties,
        couplingDicts
    );

    forAll(couplingDicts, i)
    {
        const dictionary& couplingDict = *couplingDicts[i];

        if (!couplingDict.found("conductionNetworkDomain"))
        {
            FatalErrorInFunction
                << "Conduction coupling entry #" << (i + 1)
                << " must specify conductionNetworkDomain explicitly."
                << exit(FatalError);
        }

        const word linkedConductionDomain
        (
            couplingDict.lookup("conductionNetworkDomain")
        );

        if (!conductionDomainsByName.found(linkedConductionDomain))
        {
            FatalErrorInFunction
                << "Conduction coupling entry #" << (i + 1)
                << " references conductionNetworkDomain='"
                << linkedConductionDomain
                << "', but no matching pre-myocardium domain is configured."
                << exit(FatalError);
        }

        ConductionSystemDomain& conductionDomain =
            *conductionDomainsByName[linkedConductionDomain];

        system.appendConductionCoupling
        (
            ElectroDomainCoupler::New
            (
                system.myocardium(),
                conductionDomain,
                couplingDict
            ).ptr()
        );
    }
}


void configureECGDomains
(
    electrophysicsSystem&       system,
    const electroStateProvider& stateProvider,
    const dictionary&           electroProperties
)
{
    system.endECGDomains();
    system.clearECGDomains();
    system.endECGCouplings();
    system.clearECGCouplings();

    DynamicList<word> ecgDomainNames;
    DynamicList<const dictionary*> ecgDomainDicts;

    appendSubDictionaries
    (
        electroProperties,
        "ecgDomains",
        ecgDomainNames,
        ecgDomainDicts
    );

    forAll(ecgDomainNames, i)
    {
        system.appendECGDomain
        (
            new ECGDomain
            (
                stateProvider,
                *ecgDomainDicts[i],
                ecgDomainNames[i]
            )
        );
    }
}

} // End namespace electrophysicsSystemBuilder
} // End namespace Foam

// ************************************************************************* //
