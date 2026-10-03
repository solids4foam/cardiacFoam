/*---------------------------------------------------------------------------*\
Application
    bathInterfacePatchTest

Description
    Affine heart-bath interface consistency check using the production
    exposed-phiE trace and extracellular face-conductivity helpers.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvMeshSubset.H"
#include "emptyPolyPatch.H"
#include "exposedPhiETrace.H"
#include "extracellularFaceConductivity.H"

using namespace Foam;

namespace
{

const scalar aH = 0.0409193930850485;
const scalar bH = 0.0163556464241969;
const scalar cH = 0.0220335193534877;
const scalar aB = 0.0230827345;

const tensor Ge
(
    aH, bH, 0,
    bH, cH, 0,
    0, 0, 0.0125905824877072
);
const tensor Gb(aB*tensor::I);

// The bath slope is chosen after argument parsing so that conormal flux
// is continuous for both the tangential and normal-only control cases.
scalar tangentSlope = 1.0;
scalar heartNormalSlope = 0.25;
scalar bathNormalSlope = 0.0;

scalar exactPhi(const vector& x)
{
    if (x.x() < 0)
    {
        return tangentSlope*x.y() + bathNormalSlope*x.x();
    }
    if (x.x() > 1)
    {
        return tangentSlope*x.y() + heartNormalSlope
             + bathNormalSlope*(x.x() - 1);
    }
    return tangentSlope*x.y() + heartNormalSlope*x.x();
}

wordList patchTypes(const fvMesh& mesh, const word& physicalType)
{
    wordList result(mesh.boundary().size());
    forAll(mesh.boundary(), patchI)
    {
        const word& type = mesh.boundary()[patchI].type();
        result[patchI] =
            polyPatch::constraintType(type) ? type : physicalType;
    }
    return result;
}

bool isHeart(const scalar x)
{
    return x > 0 && x < 1;
}

}


int main(int argc, char *argv[])
{
    argList::addOption("tangentSlope", "scalar", "common tangential slope");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    tangentSlope = args.getOrDefault<scalar>("tangentSlope", 1.0);
    bathNormalSlope =
        (aH*heartNormalSlope + bH*tangentSlope)/aB;

    const label zoneI = mesh.cellZones().findZoneID("myocardium");
    const label exposedPatchI = mesh.boundaryMesh().findPatchID("xMin");
    if (zoneI < 0 || exposedPatchI < 0)
    {
        FatalErrorInFunction
            << "Expected myocardium cellZone and xMin physical patch"
            << exit(FatalError);
    }

    fvMeshSubset subset(mesh);
    subset.setCellSubset(mesh.cellZones()[zoneI], exposedPatchI);
    const fvMesh& heart = subset.subMesh();
    const labelList& cellMap = subset.cellMap();
    const labelList& faceMap = subset.faceMap();

    volScalarField phiE
    (
        IOobject("affinePhiE", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("zero", dimless, 0),
        patchTypes(mesh, "fixedValue")
    );
    forAll(phiE, cellI)
    {
        phiE[cellI] = exactPhi(mesh.C()[cellI]);
    }
    forAll(phiE.boundaryField(), patchI)
    {
        if (mesh.boundary()[patchI].coupled()
         || isA<emptyPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            continue;
        }
        scalarField& values = phiE.boundaryFieldRef()[patchI];
        forAll(values, faceI)
        {
            values[faceI] = exactPhi(mesh.boundary()[patchI].Cf()[faceI]);
        }
    }
    phiE.correctBoundaryConditions();

    volTensorField sigma
    (
        IOobject("sigmaAffine", runTime.timeName(), mesh),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero),
        patchTypes(mesh, "zeroGradient")
    );
    forAll(sigma, cellI)
    {
        sigma[cellI] = isHeart(mesh.C()[cellI].x()) ? Ge : Gb;
    }
    sigma.correctBoundaryConditions();

    volTensorField heartGe
    (
        IOobject("heartGeAffine", runTime.timeName(), heart),
        heart,
        dimensionedTensor("Ge", dimless, Ge),
        patchTypes(heart, "zeroGradient")
    );
    heartGe.correctBoundaryConditions();

    volScalarField heartPhiE
    (
        IOobject("heartPhiEAffine", runTime.timeName(), heart),
        heart,
        dimensionedScalar("zero", dimless, 0),
        patchTypes(heart, "zeroGradient")
    );
    const label nExposed = useExposedPhiETrace(heartPhiE, faceMap, mesh);
    setExposedPhiETrace
    (
        heartPhiE, phiE, cellMap, faceMap, heartGe, sigma
    );

    scalar traceMax = 0;
    scalar traceSumSq = 0;
    label traceCount = 0;
    forAll(heart.boundary(), patchI)
    {
        const fvPatch& patch = heart.boundary()[patchI];
        if (patch.coupled())
        {
            continue;
        }
        forAll(patch, faceI)
        {
            const label baseFaceI = faceMap[patch.start() + faceI];
            if (baseFaceI >= mesh.nInternalFaces())
            {
                continue;
            }
            const vector& faceC = patch.Cf()[faceI];
            if (faceC.y() < 0.25 || faceC.y() > 0.75)
            {
                continue;
            }
            const scalar error =
                heartPhiE.boundaryField()[patchI][faceI]
              - exactPhi(faceC);
            traceMax = max(traceMax, mag(error));
            traceSumSq += sqr(error);
            ++traceCount;
        }
    }

    surfaceTensorField faceSigma
    (
        IOobject("faceSigmaAffine", runTime.timeName(), mesh),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );
    scalar fluxMax = 0;
    scalar fluxSumSq = 0;
    label fluxCount = 0;
    const scalar exactFlux = aH*heartNormalSlope + bH*tangentSlope;
    forAll(mesh.owner(), faceI)
    {
        const label ownerI = mesh.owner()[faceI];
        const label neighbourI = mesh.neighbour()[faceI];
        const bool ownerHeart = isHeart(mesh.C()[ownerI].x());
        const bool neighbourHeart = isHeart(mesh.C()[neighbourI].x());
        if (ownerHeart == neighbourHeart)
        {
            faceSigma[faceI] = ownerHeart ? Ge : Gb;
            continue;
        }
        const vector& faceC = mesh.Cf()[faceI];
        const scalar dOwner = mag(faceC.x() - mesh.C()[ownerI].x());
        const scalar dNeighbour = mag(mesh.C()[neighbourI].x() - faceC.x());
        faceSigma[faceI] =
            extracellularFaceConductivity::distanceWeightedHarmonic
            (
                ownerHeart ? Ge : Gb,
                neighbourHeart ? Ge : Gb,
                dOwner,
                dNeighbour
            );
        if (faceC.y() < 0.25 || faceC.y() > 0.75)
        {
            continue;
        }
        const scalar dx =
            mesh.C()[neighbourI].x() - mesh.C()[ownerI].x();
        const scalar normalSlope =
            (phiE[neighbourI] - phiE[ownerI])/dx;
        const scalar flux =
            faceSigma[faceI].xx()*normalSlope
          + faceSigma[faceI].xy()*tangentSlope;
        const scalar error = flux - exactFlux;
        fluxMax = max(fluxMax, mag(error));
        fluxSumSq += sqr(error);
        ++fluxCount;
    }
    forAll(mesh.boundary(), patchI)
    {
        faceSigma.boundaryFieldRef()[patchI] =
            sigma.boundaryField()[patchI].patchInternalField();
    }

    const volScalarField residual(fvc::laplacian(faceSigma, phiE));
    scalar stripMax = 0;
    scalar stripSumSq = 0;
    label stripCount = 0;
    forAll(residual, cellI)
    {
        const vector& c = mesh.C()[cellI];
        if (c.y() < 0.25 || c.y() > 0.75)
        {
            continue;
        }
        const scalar h = 1.0/Foam::sqrt(scalar(mesh.nCells()/3));
        if (min(mag(c.x()), mag(c.x() - 1)) > 1.5*h)
        {
            continue;
        }
        stripMax = max(stripMax, mag(residual[cellI]));
        stripSumSq += sqr(residual[cellI]);
        ++stripCount;
    }

    if (!traceCount || !fluxCount || !stripCount)
    {
        FatalErrorInFunction
            << "No interior-y interface samples found"
            << exit(FatalError);
    }

    Info<< "nExposed=" << nExposed
        << " traceCount=" << traceCount
        << " traceRms=" << Foam::sqrt(traceSumSq/traceCount)
        << " traceMax=" << traceMax
        << " fluxCount=" << fluxCount
        << " fluxRms=" << Foam::sqrt(fluxSumSq/fluxCount)
        << " fluxMax=" << fluxMax
        << " exactFlux=" << exactFlux
        << " stripCount=" << stripCount
        << " stripResidualRms=" << Foam::sqrt(stripSumSq/stripCount)
        << " stripResidualMax=" << stripMax << nl;

    return 0;
}
