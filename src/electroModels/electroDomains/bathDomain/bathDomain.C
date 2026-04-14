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

/*---------------------------------------------------------------------------*\
  NOT_COMPILED
  This file is not compiled or linked. It serves as a stub for potential
  future bath coupling development.
\*---------------------------------------------------------------------------*/

#include "bathDomain.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{

namespace
{

const dimensionSet conductivityDim
(
    pow3(dimTime) * sqr(dimCurrent)/(dimMass*dimVolume)
);

autoPtr<fvMeshSubset> createBathMeshSubset
(
    const fvMesh& supportMesh,
    const dictionary& dict
)
{
    if (dict.found("region"))
    {
        return autoPtr<fvMeshSubset>(nullptr);
    }

    if (!dict.found("cellZone"))
    {
        return autoPtr<fvMeshSubset>(nullptr);
    }

    const word cellZoneName(dict.lookup("cellZone"));
    const label zoneId = supportMesh.cellZones().findZoneID(cellZoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << "Cannot find bath cellZone '" << cellZoneName
            << "' on mesh '" << supportMesh.name() << "'."
            << exit(FatalError);
    }

    autoPtr<fvMeshSubset> subsetPtr(new fvMeshSubset(supportMesh));
    subsetPtr->setCellSubset(supportMesh.cellZones()[zoneId]);
    return subsetPtr;
}


autoPtr<fvMesh> createBathRegionMesh
(
    const fvMesh& supportMesh,
    const dictionary& dict
)
{
    if (!dict.found("region"))
    {
        return autoPtr<fvMesh>(nullptr);
    }

    const word regionName(dict.lookup("region"));

    return autoPtr<fvMesh>
    (
        new fvMesh
        (
            IOobject
            (
                regionName,
                supportMesh.time().timeName(),
                supportMesh.time(),
                IOobject::MUST_READ
            )
        )
    );
}


const fvMesh& resolveBathMesh
(
    const fvMesh& supportMesh,
    const autoPtr<fvMesh>& regionMeshPtr,
    const autoPtr<fvMeshSubset>& subsetPtr
)
{
    return
    (
        regionMeshPtr.valid()
      ? regionMeshPtr()
      :
        subsetPtr.valid() && subsetPtr->hasSubMesh()
      ? subsetPtr->subMesh()
      : supportMesh
    );
}


word resolveInterfacePatchName
(
    const fvMesh& mesh,
    const dictionary& dict,
    const autoPtr<fvMeshSubset>& subsetPtr
)
{
    const word defaultPatchName =
        (
            !dict.found("region")
         && subsetPtr.valid() && subsetPtr->hasSubMesh()
          ? fvMeshSubset::exposedPatchName
          : word::null
        );

    const word patchName =
        dict.lookupOrDefault<word>("interfacePatch", defaultPatchName);

    if (patchName.empty())
    {
        return patchName;
    }

    if (mesh.boundaryMesh().findPatchID(patchName) < 0)
    {
        FatalErrorInFunction
            << "Cannot find bath interface patch '" << patchName
            << "' on mesh '" << mesh.name() << "'."
            << exit(FatalError);
    }

    return patchName;
}


wordList bathPotentialPatchTypes
(
    const fvMesh& mesh,
    const word& interfacePatchName
)
{
    wordList patchTypes
    (
        mesh.boundary().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    if (!interfacePatchName.empty())
    {
        const label patchId = mesh.boundaryMesh().findPatchID(interfacePatchName);
        patchTypes[patchId] = fixedValueFvPatchScalarField::typeName;
    }

    return patchTypes;
}

}


defineTypeNameAndDebug(BathDomain, 0);


BathDomain::BathDomain
(
    const electroStateProvider& stateProvider,
    const fvMesh& supportMesh,
    const dictionary& dict,
    const word& domainName
)
:
    solverPtr_(BathECGSolver::New(dict)),
    stateProvider_(stateProvider),
    ownedMeshPtr_(createBathRegionMesh(supportMesh, dict)),
    meshSubsetPtr_(createBathMeshSubset(supportMesh, dict)),
    supportMesh_(supportMesh),
    mesh_(resolveBathMesh(supportMesh, ownedMeshPtr_, meshSubsetPtr_)),
    interfacePatchName_(resolveInterfacePatchName(mesh_, dict, meshSubsetPtr_)),
    interfacePatchIndex_
    (
        interfacePatchName_.empty()
      ? -1
      : mesh_.boundaryMesh().findPatchID(interfacePatchName_)
    ),
    phiE_
    (
        IOobject
        (
            domainName + "_phiE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiE", dimVoltage, 0.0),
        bathPotentialPatchTypes(mesh_, interfacePatchName_)
    ),
    interfaceSource_
    (
        IOobject
        (
            domainName + "_interfaceSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimCurrent/dimVolume, 0.0),
        "zeroGradient"
    ),
    sigmaBath_
    (
        IOobject
        (
            domainName + "_sigmaBath",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("sigmaBath", conductivityDim, 0.0),
        "zeroGradient"
    )
{
    if (!sigmaBath_.headerOk())
    {
        const scalar sigmaBathValue =
            dict.lookupOrDefault<scalar>("sigmaBath", -1.0);

        if (sigmaBathValue < 0.0)
        {
            FatalErrorInFunction
                << "BathDomain '" << domainName
                << "' requires either a field file named '"
                << sigmaBath_.name() << "' or a scalar 'sigmaBath' entry in "
                << dict.dictName() << "."
                << exit(FatalError);
        }

        sigmaBath_ =
            dimensionedScalar("sigmaBath", conductivityDim, sigmaBathValue);
    }

    Info<< "Constructed BathDomain '" << domainName
        << "' on mesh '" << mesh_.name() << "'";

    if (ownedMeshPtr_.valid())
    {
        Info<< " from region '"
            << dict.lookupOrDefault<word>("region", word::null) << "'";
    }
    else if (usesSubsetMesh())
    {
        Info<< " from cellZone '"
            << dict.lookupOrDefault<word>("cellZone", word::null) << "'";
    }

    if (interfacePatchIndex_ >= 0)
    {
        Info<< " with interface patch '" << interfacePatchName_ << "'";
    }

    Info<< "." << nl << endl;

    Info<< "BathDomain '" << domainName << "' boundary patches:" << nl;
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];
        Info<< "  [" << patchi << "] " << patch.name()
            << " : " << patch.size() << " faces" << nl;
    }

    if (interfacePatchIndex_ >= 0)
    {
        const polyPatch& interfacePatch = mesh_.boundaryMesh()[interfacePatchIndex_];

        Info<< "BathDomain '" << domainName << "' selected interface patch '"
            << interfacePatch.name() << "' has "
            << interfacePatch.size() << " faces." << nl;

        if (!interfacePatch.size())
        {
            WarningInFunction
                << "Selected interface patch '" << interfacePatch.name()
                << "' is empty on bath mesh '" << mesh_.name() << "'."
                << nl;
        }
    }

    Info<< endl;
}


void BathDomain::advance
(
    scalar t0,
    scalar dt
)
{
    solverPtr_->solve(*this, t0, dt);
}

} // End namespace Foam

// ************************************************************************* //
