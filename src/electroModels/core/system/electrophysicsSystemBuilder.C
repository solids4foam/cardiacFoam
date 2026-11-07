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
#include "heartPhiETrace.H"
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
    const word&              solverType,
    const dictionary&        electroProperties,
    scalar                   initialDeltaT,
    const torsoCoupling::type coupling
)
{
    system.setMyocardium
    (
        myocardiumDomainInterface::New
        (
            mesh,
            solverType,
            electroProperties,
            initialDeltaT,
            coupling == torsoCoupling::type::bath
        ).ptr()
    );

    // The heart phiE of a bidomain at the heart-torso wall:
    // bathHeartPhiETrace is required with a torso and rejected without one
    if (system.myocardium().hasExtracellularPotential())
    {
        (void)heartPhiETrace::read(electroProperties);
    }
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
    const dictionary&     electroProperties,
    const torsoCoupling::type coupling
)
{
    system.clearPotentialDomain();

    if (coupling != torsoCoupling::type::bath)
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

    if (!system.myocardium().hasExtracellularPotential())
    {
        FatalIOErrorInFunction(electroProperties)
            << "torsoCoupling bath solves phiE of the heart and the torso "
            << "together, which needs myocardiumSolver bidomainSolver; this "
            << "myocardium has no extracellular potential. Use torsoCoupling "
            << "leadField for a one-way torso ECG of a monodomainSolver heart."
            << exit(FatalIOError);
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
    const myocardiumDomainInterface& myocardium,
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

            ecgDomainNames.append(e.keyword());
            ecgDomainDicts.append(&e.dict());
        }
    }

    HashTable<ecgDomain*> ecgDomainsByName(ecgDomainNames.size());

    // The output file (pseudoECG.dat, torsoECG.dat, eikonalECG.dat) is
    // named after the ecgSolver, not the domain: two domains with the same
    // ecgSolver write the same file
    HashTable<word> ecgDomainBySolverType(ecgDomainNames.size());

    forAll(ecgDomainNames, i)
    {
        const word& domainName = ecgDomainNames[i];
        const dictionary& domainDict = *ecgDomainDicts[i];

        // The solver says what it reads and which myocardium it needs
        autoPtr<ecgSolver> solverPtr(ecgSolver::New(domainDict));
        const word ecgSolverType(solverPtr->type());

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

        solverPtr->checkMyocardium(myocardium, domainName);

        const electroStateProvider* stateProviderPtr =
            solverPtr->source() == ecgSolver::stateSource::bath
          ? potentialStateProviderPtr
          : &myocardium;

        if (!stateProviderPtr)
        {
            FatalErrorInFunction
                << "ECG domain '" << domainName << "' selects ecgSolver "
                << ecgSolverType << ", which reads the bath potential "
                << "(torsoCoupling bath), but no bath domain is configured."
                << exit(FatalError);
        }

        ecgDomain* domainPtr =
            new ecgDomain
            (
                *stateProviderPtr,
                std::move(solverPtr),
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
