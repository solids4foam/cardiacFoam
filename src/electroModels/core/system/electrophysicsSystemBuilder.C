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
#include "torsoCoupling.H"
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

    if (torsoCoupling::read(electroProperties) != torsoCoupling::type::bath)
    {
        return;
    }

    const word solverType(electroProperties.get<word>("myocardiumSolver"));

    if (solverType != "bidomainSolver")
    {
        FatalIOErrorInFunction(electroProperties)
            << "torsoCoupling bath solves phiE of the heart and the torso "
            << "together, which needs myocardiumSolver bidomainSolver; "
            << "myocardiumSolver is '" << solverType << "'. Use torsoCoupling "
            << "leadField for a one-way torso ECG of a monodomainSolver heart."
            << exit(FatalIOError);
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
            electroProperties
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

    // One network: its root stimulus is the start of activation, and its
    // outputs (purkinjeNetwork.dat, ...) have fixed names
    if (conductionDomainNames.size() > 1)
    {
        FatalIOErrorInFunction
        (
            electroProperties.subDict("conductionNetworkDomains")
        )
            << "conductionNetworkDomains holds "
            << conductionDomainNames.size() << " domains "
            << conductionDomainNames
            << ", but only one conduction network domain is supported."
            << exit(FatalIOError);
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
                domainName,
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

    // A network that no coupling links to the myocardium would be advanced
    // without ever activating it
    if (couplingDicts.empty())
    {
        FatalIOErrorInFunction(electroProperties)
            << "conductionNetworkDomains holds " << conductionDomainNames
            << ", but domainCouplings holds no coupling that links it to "
            << "the myocardium." << exit(FatalIOError);
    }

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
    system.endECGDomains();
    system.clearECGDomains();

    DynamicList<word> ecgDomainNames;
    DynamicList<const dictionary*> ecgDomainDicts;
    const dictionary* sharedElectrodePositionsPtr = nullptr;

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

            // An ECG reads the myocardium (and bath) state directly, and
            // torsoCoupling couples the torso, so there is nothing to couple
            // an ECG domain through
            if (e.dict().found("coupling"))
            {
                FatalIOErrorInFunction(e.dict())
                    << "ECG domain '" << e.keyword()
                    << "' has a coupling entry, but ECG domains take no "
                    << "coupling." << exit(FatalIOError);
            }

            ecgDomainNames.append(e.keyword());
            ecgDomainDicts.append(&e.dict());
        }
    }

    HashTable<ecgDomain*> ecgDomainsByName(ecgDomainNames.size());

    // The output file (pseudoECG.dat, torsoECG.dat, eikonalECG.dat) is
    // named after the ecgSolver, not the domain: two domains with the same
    // ecgSolver write the same file
    HashTable<word> ecgDomainBySolverType(ecgDomainNames.size());

    const word myocardiumSolver
    (
        electroProperties.get<word>("myocardiumSolver")
    );
    const torsoCoupling::type coupling = torsoCoupling::read(electroProperties);

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

        if (ecgDomainBySolverType.found(ecgSolverType))
        {
            WarningInFunction
                << "ECG domains '" << ecgDomainBySolverType[ecgSolverType]
                << "' and '" << domainName << "' both select ecgSolver "
                << ecgSolverType << " and write the same output file in "
                << "postProcessing/; each write overwrites the other's."
                << endl;
        }
        else
        {
            ecgDomainBySolverType.insert(ecgSolverType, domainName);
        }

        const electroStateProvider* stateProviderPtr = nullptr;

        if (ecgSolverType == "pseudoECG")
        {
            stateProviderPtr = &myocardiumStateProvider;
        }
        else if (ecgSolverType == "eikonalECG")
        {
            if (myocardiumSolver != "eikonalSolver")
            {
                FatalErrorInFunction
                    << "ECG domain '" << domainName
                    << "' selects ecgSolver eikonalECG, but "
                    << "myocardiumSolver is '" << myocardiumSolver
                    << "'. eikonalECG requires myocardiumSolver "
                    << "eikonalSolver because it reconstructs Vm from "
                    << "the activation-time field."
                    << exit(FatalError);
            }

            stateProviderPtr = &myocardiumStateProvider;
        }
        else if (ecgSolverType == "torsoECG")
        {
            if (coupling == torsoCoupling::type::none)
            {
                FatalIOErrorInFunction(domainDict)
                    << "ECG domain '" << domainName
                    << "' selects ecgSolver torsoECG, but electroProperties "
                    << "has no torso. Add the region outside the heart and "
                    << "how its potential is obtained:" << nl
                    << "torso" << nl
                    << "{" << nl
                    << "    cellZones         (...);" << nl
                    << "    conductivityField <name>;" << nl
                    << "    ..." << nl
                    << "}" << nl
                    << "torsoCoupling bath | leadField;"
                    << exit(FatalIOError);
            }

            // bath: the bath domain's phiE; leadField: the myocardium's Vm
            stateProviderPtr =
                coupling == torsoCoupling::type::bath
              ? potentialStateProviderPtr
              : &myocardiumStateProvider;

            if (!stateProviderPtr)
            {
                FatalErrorInFunction
                    << "ECG domain '" << domainName << "' selects ecgSolver "
                    << "torsoECG with torsoCoupling bath, but no bath domain "
                    << "is configured."
                    << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "ECG domain '" << domainName
                << "' selects ecgSolver '" << ecgSolverType
                << "', but provider routing is only defined for "
                << "pseudoECG, eikonalECG and torsoECG."
                << exit(FatalError);
        }

        ecgDomain* domainPtr =
            new ecgDomain
            (
                *stateProviderPtr,
                domainDict,
                domainName,
                sharedElectrodePositionsPtr
            );

        ecgDomainsByName.insert(domainName, domainPtr);
        system.appendECGDomain(domainPtr);
    }
}

} // End namespace electrophysicsSystemBuilder
} // End namespace Foam

// ************************************************************************* //
