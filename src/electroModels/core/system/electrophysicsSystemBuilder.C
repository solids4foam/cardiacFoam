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
#include "extracellularPotentialDomain.H"
#include "electrophysicsAdvanceScheme.H"
#include "error.H"

#include "DynamicList.H"
#include "HashTable.H"

namespace Foam
{
namespace electrophysicsSystemBuilder
{

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


word myocardiumSolverType(const dictionary& electroProperties)
{
    if (electroProperties.found("myocardiumSolver"))
    {
        return word(electroProperties.lookup("myocardiumSolver"));
    }

    const word coeffDictName = electroProperties.dictName();

    return coeffDictName.endsWith("Coeffs")
      ? word(coeffDictName.substr(0, coeffDictName.size() - 6))
      : word("unset");
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


void configureBathPotentialDomain
(
    electrophysicsSystem& system,
    const fvMesh&         mesh,
    const dictionary&     electroProperties
)
{
    system.clearPotentialDomain();

    if (!electroProperties.found("bathPotentialDomain"))
    {
        return;
    }

    if (!system.hasMyocardium())
    {
        FatalErrorInFunction
            << "configureBathPotentialDomain requires the myocardium domain "
            << "to be configured first."
            << exit(FatalError);
    }

    system.setPotentialDomain
    (
        new extracellularPotentialDomain
        (
            mesh,
            system.myocardium(),
            electroProperties.subDict("bathPotentialDomain")
        )
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

    HashTable<conductionSystemDomain*> conductionDomainsByName
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

        autoPtr<conductionSystemDomain> conductionDomain
        (
            conductionSystemDomain::New
            (
                mesh,
                *conductionDomainDicts[i],
                initialDeltaT
            )
        );

        conductionSystemDomain* domainPtr = conductionDomain.ptr();
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

        conductionSystemDomain& conductionDomain =
            *conductionDomainsByName[linkedConductionDomain];

        system.appendConductionCoupling
        (
            electroDomainCoupler::New
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
    const electroStateProvider& myocardiumStateProvider,
    const electroStateProvider* potentialStateProviderPtr,
    const dictionary&           electroProperties
)
{
    system.endECGCouplings();
    system.clearECGCouplings();
    system.endECGDomains();
    system.clearECGDomains();

    DynamicList<word> ecgDomainNames;
    DynamicList<const dictionary*> ecgDomainDicts;
    const dictionary* sharedElectrodePositionsPtr = nullptr;
    const dictionary* manufacturedBidomainPtr =
        electroProperties.findDict("manufacturedBidomain");
    const dictionary* bathPotentialDomainPtr =
        electroProperties.findDict("bathPotentialDomain");

    if (electroProperties.found("ecgDomains"))
    {
        const dictionary& ecgDomainsDict =
            electroProperties.subDict("ecgDomains");

        sharedElectrodePositionsPtr =
            ecgDomainsDict.findDict("electrodePositions");

        forAllConstIter(dictionary, ecgDomainsDict, iter)
        {
            const entry& e = iter();

            if (!e.isDict() || e.keyword() == "electrodePositions")
            {
                continue;
            }

            ecgDomainNames.append(e.keyword());
            ecgDomainDicts.append(&e.dict());
        }
    }

    HashTable<ecgDomain*> ecgDomainsByName(ecgDomainNames.size());

    forAll(ecgDomainNames, i)
    {
        const word& domainName = ecgDomainNames[i];
        const dictionary& domainDict = *ecgDomainDicts[i];
        const word ecgSolverType
        (
            domainDict.lookupOrDefault<word>("ecgSolver", "pseudoECG")
        );

        if (ecgDomainsByName.found(domainName))
        {
            FatalErrorInFunction
                << "Duplicate ECG domain name '" << domainName
                << "' while configuring post-myocardium domains."
                << exit(FatalError);
        }

        const electroStateProvider* stateProviderPtr = nullptr;

        if (ecgSolverType == "pseudoECG")
        {
            stateProviderPtr = &myocardiumStateProvider;
        }
        else if (ecgSolverType == "torsoECG")
        {
            const word solverType(myocardiumSolverType(electroProperties));

            if (solverType != "bidomainSolver")
            {
                FatalErrorInFunction
                    << "ECG domain '" << domainName
                    << "' selects ecgSolver torsoECG, but "
                    << "myocardiumSolver is '" << solverType
                    << "'. torsoECG requires myocardiumSolver "
                    << "bidomainSolver because it samples the "
                    << "extracellular potential phiE."
                    << exit(FatalError);
            }

            if (!potentialStateProviderPtr)
            {
                FatalErrorInFunction
                    << "ECG domain '" << domainName
                    << "' selects ecgSolver torsoECG, but no "
                    << "bathPotentialDomain is configured inside "
                    << "bidomainSolverCoeffs. Add:" << nl
                    << "bidomainSolverCoeffs" << nl
                    << "{" << nl
                    << "    bathPotentialDomain" << nl
                    << "    {" << nl
                    << "        bathCellZones (...);" << nl
                    << "        ..." << nl
                    << "    }" << nl
                    << "}" << exit(FatalError);
            }

            stateProviderPtr = potentialStateProviderPtr;
        }
        else
        {
            FatalErrorInFunction
                << "ECG domain '" << domainName
                << "' selects ecgSolver '" << ecgSolverType
                << "', but provider routing is only defined for "
                << "pseudoECG and torsoECG."
                << exit(FatalError);
        }

        ecgDomain* domainPtr =
            new ecgDomain
            (
                *stateProviderPtr,
                domainDict,
                domainName,
                sharedElectrodePositionsPtr,
                manufacturedBidomainPtr,
                bathPotentialDomainPtr
            );

        ecgDomainsByName.insert(domainName, domainPtr);
        system.appendECGDomain(domainPtr);
    }

    forAll(ecgDomainNames, i)
    {
        const dictionary& domainDict = *ecgDomainDicts[i];

        if (!domainDict.found("coupling"))
        {
            continue;
        }

        if (!system.hasMyocardium())
        {
            FatalErrorInFunction
                << "ECG domain '" << ecgDomainNames[i]
                << "' configures a coupling block, but no myocardium domain "
                << "is available as the primary coupling endpoint."
                << exit(FatalError);
        }

        system.appendECGCoupling
        (
            electroDomainCoupler::New
            (
                system.myocardium(),
                *ecgDomainsByName[ecgDomainNames[i]],
                domainDict.subDict("coupling")
            ).ptr()
        );
    }
}

} // End namespace electrophysicsSystemBuilder
} // End namespace Foam

// ************************************************************************* //
