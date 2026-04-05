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
#include "myocardiumDomain.H"
#include "electroDomainCoupler.H"
#include "electrophysicsAdvanceScheme.H"
#include "error.H"
#include "ionicModel.H"
#include "ecgDomain.H"
#include "myocardiumDomain.H"
#include "myocardiumSolver.H"
#include "monodomainSolver.H"

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
    const word&       entryName,
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


void collectConductionSystemDomainDicts
(
    const dictionary& electroProperties,
    DynamicList<word>& names,
    DynamicList<const dictionary*>& dicts
)
{
    if
    (
        appendSubDictionaries
        (
            electroProperties,
            "conductionNetworkDomains",
            names,
            dicts
        )
    )
    {
        return;
    }

    if (electroProperties.found("auxiliaryDomains"))
    {
        FatalErrorInFunction
            << "Unsupported legacy key 'auxiliaryDomains'. "
            << "Use 'conductionNetworkDomains' instead."
            << exit(FatalError);
    }

    if (electroProperties.found("purkinjeNetwork"))
    {
        FatalErrorInFunction
            << "Unsupported legacy top-level key 'purkinjeNetwork'. "
            << "Use 'conductionNetworkDomains { <name> { ... } }' instead."
            << exit(FatalError);
    }
}


void collectConductionCouplingDicts
(
    const dictionary& electroProperties,
    DynamicList<word>& names,
    DynamicList<const dictionary*>& dicts
)
{
    if
    (
        appendSubDictionaries
        (
            electroProperties,
            "domainCouplings",
            names,
            dicts
        )
    )
    {
        return;
    }

    if (electroProperties.found("purkinjeNetwork"))
    {
        FatalErrorInFunction
            << "Unsupported legacy coupling configuration under "
            << "'purkinjeNetwork'. Use 'domainCouplings { <name> { ... } }'."
            << exit(FatalError);
    }
}

} // End anonymous namespace


void configureMyocardiumDomain
(
    electrophysicsSystem&   system,
    const fvMesh&           mesh,
    const dictionary&       electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList&          postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    ionicModel&              ionicModel,
    electroVerificationModel* verificationModelPtr,
    scalar                  initialDeltaT
)
{
    // Derive the solver type from the Coeffs subdict name
    // ("monodomainSolverCoeffs" → "monodomainSolver").
    const word& cn = electroProperties.dictName();
    const word solverType =
        cn.endsWith("Coeffs") ? word(cn.substr(0, cn.size() - 6)) : cn;

    system.setMyocardium
    (
        new MyocardiumDomain
        (
            mesh,
            electroProperties,
            outFields,
            postProcessFieldNames,
            postProcessFields,
            ionicModel,
            verificationModelPtr,
            myocardiumSolver::New(mesh, solverType, electroProperties)
        )
    );

    configureAdvanceScheme(system, electroProperties);
    configureConductionSystemDomain
    (
        system, mesh, electroProperties, initialDeltaT
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


void configureConductionSystemDomain
(
    electrophysicsSystem& system,
    const fvMesh&         mesh,
    const dictionary&     electroProperties,
    scalar                initialDeltaT
)
{
    const ConductionSystemDomainContext context(mesh.time(), &mesh);

    system.clearUpstreamDomain();
    system.clearUpstreamCouplingModel();

    DynamicList<word> networkNames;
    DynamicList<const dictionary*> networkDicts;
    collectConductionSystemDomainDicts
    (
        electroProperties,
        networkNames,
        networkDicts
    );

    if (networkNames.empty())
    {
        return;
    }

    HashTable<ConductionSystemDomain*> networkDomainsByName
    (
        networkNames.size()
    );

    forAll(networkNames, i)
    {
        const word& networkName = networkNames[i];

        if (networkDomainsByName.found(networkName))
        {
            FatalErrorInFunction
                << "Duplicate conduction domain name '" << networkName
                << "' while configuring pre-primary domains."
                << exit(FatalError);
        }

        autoPtr<ConductionSystemDomain> networkDomain =
            ConductionSystemDomain::New
            (
                context,
                *networkDicts[i],
                initialDeltaT
            );

        ConductionSystemDomain* networkPtr = networkDomain.ptr();
        networkDomainsByName.insert(networkName, networkPtr);
        system.appendUpstreamDomain(networkPtr);
    }

    DynamicList<word> couplingNames;
    DynamicList<const dictionary*> couplingDicts;
    collectConductionCouplingDicts
    (
        electroProperties,
        couplingNames,
        couplingDicts
    );

    forAll(couplingNames, i)
    {
        const dictionary& couplingDict = *couplingDicts[i];

        const word primaryName =
            couplingDict.lookupOrDefault<word>("primaryDomain", "myocardium");

        if (primaryName != "myocardium")
        {
            FatalErrorInFunction
                << "Only primaryDomain=myocardium is currently supported, got '"
                << primaryName << "' in coupling entry '" << couplingNames[i]
                << "'."
                << exit(FatalError);
        }

        word linkedNetwork =
            couplingDict.lookupOrDefault<word>
            (
                "conductionNetworkDomain",
                word::null
            );

        if (linkedNetwork.empty())
        {
            if (networkNames.size() == 1)
            {
                linkedNetwork = networkNames[0];
            }
            else
            {
                FatalErrorInFunction
                    << "Coupling entry '" << couplingNames[i]
                    << "' must specify conductionNetworkDomain when multiple "
                    << "pre-primary domains are configured."
                    << exit(FatalError);
            }
        }

        if (!networkDomainsByName.found(linkedNetwork))
        {
            FatalErrorInFunction
                << "Coupling entry '" << couplingNames[i]
                << "' references conductionNetworkDomain='"
                << linkedNetwork << "', but no matching pre-primary domain is "
                << "configured."
                << exit(FatalError);
        }

        ConductionSystemDomain& networkDomain =
            *networkDomainsByName[linkedNetwork];

        system.appendUpstreamCouplingModel
        (
            ElectroDomainCoupler::New
            (
                system.myocardium(),
                networkDomain,
                couplingDict
            ).ptr()
        );
    }
}


void configureECGDomain
(
    electrophysicsSystem&       system,
    const electroStateProvider& stateProvider,
    const dictionary&           electroProperties
)
{
    system.endDownstreamDomain();
    system.clearDownstreamDomain();
    system.endDownstreamCouplingModel();
    system.clearDownstreamCouplingModel();

    DynamicList<word> ecgDomainNames;
    DynamicList<const dictionary*> ecgDomainDicts;

    if
    (
        appendSubDictionaries
        (
            electroProperties,
            "ecgDomains",
            ecgDomainNames,
            ecgDomainDicts
        )
    )
    {
        // configured via ecgDomains
    }
    else if (electroProperties.found("ECG"))
    {
        FatalErrorInFunction
            << "Unsupported legacy key 'ECG'. "
            << "Use 'ecgDomains { <name> { ... } }' instead."
            << exit(FatalError);
    }

    forAll(ecgDomainNames, i)
    {
        system.appendDownstreamDomain
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
