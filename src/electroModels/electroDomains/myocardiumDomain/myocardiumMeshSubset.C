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

#include "myocardiumMeshSubset.H"
#include "emptyPolyPatch.H"
#include "error.H"

namespace Foam
{

label myocardiumExposedFacesPatch
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties
)
{
    const polyBoundaryMesh& patches = supportMesh.boundaryMesh();
    label patchI = -1;

    if (electroProperties.found("exposedFacesPatch"))
    {
        const word patchName(electroProperties.get<word>("exposedFacesPatch"));
        patchI = patches.findPatchID(patchName);

        if (patchI < 0)
        {
            FatalIOErrorInFunction(electroProperties)
                << "exposedFacesPatch '" << patchName << "' is not a patch of "
                << "mesh '" << supportMesh.name() << "'."
                << exit(FatalIOError);
        }
    }
    else
    {
        forAll(patches, candidateI)
        {
            const polyPatch& pp = patches[candidateI];
            if (!isA<emptyPolyPatch>(pp) && !pp.coupled())
            {
                patchI = candidateI;
                break;
            }
        }

        if (patchI < 0)
        {
            FatalIOErrorInFunction(electroProperties)
                << "Mesh '" << supportMesh.name() << "' has no patch that is "
                << "neither empty nor coupled to receive the faces exposed by "
                << "the myocardium cellZone. Name one with exposedFacesPatch."
                << exit(FatalIOError);
        }
    }

    const polyPatch& pp = patches[patchI];
    if (isA<emptyPolyPatch>(pp) || pp.coupled())
    {
        FatalIOErrorInFunction(electroProperties)
            << "exposedFacesPatch '" << pp.name() << "' is " << pp.type()
            << "; the exposed myocardium faces need a patch that is neither "
            << "empty nor coupled."
            << exit(FatalIOError);
    }

    Info<< "Myocardium cellZone: exposed faces go into patch '" << pp.name()
        << "' (" << pp.type() << ")." << endl;

    return patchI;
}


autoPtr<fvMeshSubset> createMyocardiumMeshSubset
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties
)
{
    if (!electroProperties.found("cellZone"))
    {
        return autoPtr<fvMeshSubset>(nullptr);
    }

    const word cellZoneName(electroProperties.lookup("cellZone"));
    const label zoneId = supportMesh.cellZones().findZoneID(cellZoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << "Cannot find myocardium cellZone '" << cellZoneName
            << "' on mesh '" << supportMesh.name() << "'."
            << exit(FatalError);
    }

    autoPtr<fvMeshSubset> subsetPtr(new fvMeshSubset(supportMesh));
    subsetPtr->setCellSubset
    (
        supportMesh.cellZones()[zoneId],
        myocardiumExposedFacesPatch(supportMesh, electroProperties)
    );
#if OPENFOAM < 2406
    // OpenFOAM before 2406 copies the base fvSchemes and fvSolution onto the
    // sub-mesh without parsing them, so read them for the scheme and solver
    // lookups to resolve.
    fvMesh& subMesh = subsetPtr->subMesh();
    subMesh.fvSchemes::readOpt(IOobject::MUST_READ);
    subMesh.fvSchemes::read();
    subMesh.fvSolution::readOpt(IOobject::MUST_READ);
    subMesh.fvSolution::read();
#endif
    return subsetPtr;
}


} // End namespace Foam

// ************************************************************************* //
