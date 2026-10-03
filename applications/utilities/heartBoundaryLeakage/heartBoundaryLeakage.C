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

Application
    heartBoundaryLeakage

Description
    Net current through the sealed heart boundary.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvMeshSubset.H"
#include "timeSelector.H"
#include "emptyPolyPatch.H"
#include "insulatedFaceConductivity.H"

using namespace Foam;

namespace
{

label defaultExposedPatch(const fvMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchI)
    {
        if (!isA<emptyPolyPatch>(patches[patchI]) && !patches[patchI].coupled())
        {
            return patchI;
        }
    }
    return -1;
}


tmp<volTensorField> readConductivity(const fvMesh& mesh, const word& name)
{
    IOobject io
    (
        name, "0", mesh, IOobject::MUST_READ, IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (io.typeHeaderOk<volSymmTensorField>(true))
    {
        const volSymmTensorField G(io, mesh);
        return tmp<volTensorField>::New
        (
            IOobject(name + "Tensor", "0", mesh, IOobject::NO_READ,
                IOobject::NO_WRITE, IOobject::NO_REGISTER),
            G & tensor(I)
        );
    }
    return tmp<volTensorField>::New(io, mesh);
}


wordList subPatchTypes(const fvMesh& sub, const word& physicalType)
{
    wordList types(sub.boundary().size());
    forAll(sub.boundary(), patchI)
    {
        const word& t = sub.boundary()[patchI].type();
        types[patchI] = polyPatch::constraintType(t) ? t : physicalType;
    }
    return types;
}


tmp<volScalarField> subPotential
(
    const fvMesh& sub,
    const labelList& cellMap,
    const volScalarField& base,
    const bool analyticXY,
    const wordList& types
)
{
    tmp<volScalarField> tpsi = tmp<volScalarField>::New
    (
        IOobject(base.name() + "Sub", sub.time().timeName(), sub),
        sub,
        dimensionedScalar(base.dimensions(), Zero),
        types
    );
    volScalarField& psi = tpsi.ref();
    const vectorField& C = sub.C().primitiveField();

    forAll(psi, cellI)
    {
        psi[cellI] =
            analyticXY
          ? C[cellI].x()*C[cellI].y()
          : base.primitiveField()[cellMap[cellI]];
    }
    return tpsi;
}


tmp<volScalarField> sealedPotential
(
    const fvMesh& sub,
    const labelList& cellMap,
    const volScalarField& base,
    const bool analyticXY
)
{
    tmp<volScalarField> tpsi = subPotential
    (
        sub, cellMap, base, analyticXY, subPatchTypes(sub, "zeroGradient")
    );
    tpsi.ref().correctBoundaryConditions();
    return tpsi;
}


tmp<volScalarField> globalTracePotential
(
    const fvMesh& sub,
    const labelList& cellMap,
    const labelList& faceMap,
    const volScalarField& base,
    const volSymmTensorField* GePtr,
    const volScalarField* sigmaBathPtr
)
{
    const fvMesh& baseMesh = base.mesh();
    wordList types(subPatchTypes(sub, "zeroGradient"));
    forAll(sub.boundary(), patchI)
    {
        const fvPatch& p = sub.boundary()[patchI];
        if (polyPatch::constraintType(p.type()))
        {
            continue;
        }
        forAll(p, i)
        {
            if (faceMap[p.start() + i] < baseMesh.nInternalFaces())
            {
                types[patchI] = "fixedValue";
            }
        }
    }

    tmp<volScalarField> tpsi =
        subPotential(sub, cellMap, base, false, types);
    volScalarField& psi = tpsi.ref();

    const scalarField& w = baseMesh.weights().primitiveField();
    forAll(psi.boundaryField(), patchI)
    {
        if (types[patchI] != "fixedValue")
        {
            continue;
        }
        const fvPatch& p = sub.boundary()[patchI];
        fvPatchScalarField& pf = psi.boundaryFieldRef()[patchI];
        forAll(pf, i)
        {
            const label f = faceMap[p.start() + i];
            if (f < baseMesh.nInternalFaces() && GePtr)
            {
                const label hc = cellMap[p.faceCells()[i]];
                const label bc =
                    baseMesh.owner()[f] == hc
                  ? baseMesh.neighbour()[f]
                  : baseMesh.owner()[f];
                const vector n = baseMesh.Sf()[f]/baseMesh.magSf()[f];
                const scalar dH =
                    mag((baseMesh.Cf()[f] - baseMesh.C()[hc]) & n);
                const scalar dB =
                    mag((baseMesh.Cf()[f] - baseMesh.C()[bc]) & n);
                const scalar kH = (n & ((*GePtr)[hc] & n))/dH;
                const scalar kB = (*sigmaBathPtr)[bc]/dB;
                pf[i] = (kH*base[hc] + kB*base[bc])/(kH + kB);
            }
            else if (f < baseMesh.nInternalFaces())
            {
                pf[i] =
                    w[f]*base[baseMesh.owner()[f]]
                  + (1 - w[f])*base[baseMesh.neighbour()[f]];
            }
            else
            {
                const label bp = baseMesh.boundaryMesh().whichPatch(f);
                pf[i] =
                    base.boundaryField()[bp]
                    [f - baseMesh.boundaryMesh()[bp].start()];
            }
        }
    }
    psi.correctBoundaryConditions();
    return tpsi;
}


tmp<volScalarField> conormalPotential
(
    const fvMesh& sub,
    const labelList& cellMap,
    const volScalarField& base,
    const bool analyticXY,
    const volTensorField& G,
    const label nIter
)
{
    tmp<volScalarField> tpsi = subPotential
    (
        sub, cellMap, base, analyticXY, subPatchTypes(sub, "conormalZeroFlux")
    );
    volScalarField& psi = tpsi.ref();
    psi.correctBoundaryConditions();

    const volTensorField Gc
    (
        IOobject("conductivity", sub.time().timeName(), sub),
        G
    );
    volVectorField gradPsi
    (
        IOobject("grad(" + psi.name() + ")", sub.time().timeName(), sub),
        sub,
        dimensionedVector(psi.dimensions()/dimLength, Zero)
    );
    gradPsi = fvc::grad(psi);
    for (label iter = 0; iter < nIter; ++iter)
    {
        psi.correctBoundaryConditions();
        gradPsi = fvc::grad(psi);
    }
    return tpsi;
}

}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption("zone", "word", "myocardium | none");
    argList::addOption("exposedPatch", "word", "patch of exposed faces");
    argList::addOption("conductivity", "word", "0/ conductivity field");
    argList::addOption("fields", "wordList", "potentials");
    argList::addOption("refPoint", "point", "phiERefPoint");
    argList::addOption("refValue", "scalar", "phiEReferenceValue");
    argList::addOption("sigmaBath", "word", "0/ bath conductivity field");
    argList::addOption("analyticField", "word", "xy");
    argList::addBoolOption("insulated", "insulatedFaceConductivity");
    argList::addBoolOption("insulationDelta", "max|laplacian sealed - laplacian| per field");
    argList::addBoolOption("gradientBias", "wall-cell grad(phiE) bias");
    argList::addOption("phiETrace", "word", "zeroGradient | global | conductivityWeighted");
    argList::addOption("conductivityExtracellular", "word", "0/ extracellular conductivity field");
    argList::addOption("exactPhiE", "word", "bathFDA: exact grad(phiE) reference for -gradientBias");
    argList::addOption("mmsK", "scalar", "bathFDA k");
    argList::addOption("mmsAlpha", "scalar", "bathFDA alpha");
    argList::addOption("mmsSe", "scalar", "bathFDA s_e");
    argList::addOption("conormalIterations", "label", "conormalZeroFlux passes");
    argList::addOption("restValue", "scalar", "Vm rest value for -restDrift");
    argList::addBoolOption("restDrift", "max|Vm - restValue| on wall and interior heart cells");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const word zoneName(args.getOrDefault<word>("zone", "myocardium"));
    const word conductivityName
    (
        args.getOrDefault<word>("conductivity", "ConductivityIntracellular")
    );
    const wordList fieldNames
    (
        args.getOrDefault<wordList>("fields", wordList({"Vm", "phiE"}))
    );
    const bool analyticXY(args.getOrDefault<word>("analyticField", "") == "xy");
    const bool insulated(args.found("insulated"));
    const label conormalIterations
    (
        args.getOrDefault<label>("conormalIterations", 0)
    );

    autoPtr<fvMeshSubset> subsetterPtr;
    label exposedI = -1;
    if (zoneName != "none")
    {
        const label zoneI = mesh.cellZones().findZoneID(zoneName);
        if (zoneI < 0)
        {
            FatalErrorInFunction
                << "cellZone " << zoneName << " not found" << exit(FatalError);
        }
        exposedI =
            args.found("exposedPatch")
          ? mesh.boundaryMesh().findPatchID(args.get<word>("exposedPatch"))
          : defaultExposedPatch(mesh);

        subsetterPtr.reset
        (
            new fvMeshSubset
            (
                mesh,
                bitSet(mesh.nCells(), mesh.cellZones()[zoneI]),
                exposedI
            )
        );
    }
    const fvMesh& sub = subsetterPtr ? subsetterPtr->subMesh() : mesh;
    labelList wholeMeshMap;
    if (!subsetterPtr)
    {
        wholeMeshMap = Foam::identity(mesh.nCells());
    }
    const labelList& cellMap =
        subsetterPtr ? subsetterPtr->cellMap() : wholeMeshMap;

    const tmp<volTensorField> tGbase(readConductivity(mesh, conductivityName));
    const volTensorField G
    (
        IOobject(conductivityName + "Sub", "0", sub, IOobject::NO_READ,
            IOobject::NO_WRITE, IOobject::NO_REGISTER),
        subsetterPtr
      ? subsetterPtr->interpolate(tGbase())()
      : tGbase()
    );

    label refCell = -1;
    scalar dRef = 0;
    const bool haveRef = args.found("refPoint");
    const scalar refValue(args.getOrDefault<scalar>("refValue", 0));

    Info<< "# time";
    forAll(fieldNames, i)
    {
        Info<< " L_" << fieldNames[i] << " Lnorm_" << fieldNames[i];
    }
    Info<< " L";
    if (args.found("restDrift"))
    {
        Info<< " driftWall driftInner";
    }
    if (args.found("gradientBias"))
    {
        Info<< " gradBiasWall gradBiasWallMax gradBiasInner";
        if (args.found("exactPhiE"))
        {
            Info<< " globalBiasWall";
        }
    }
    if (haveRef)
    {
        Info<< " d_r phiRefPredicted phiRefSolved relErr";
    }
    Info<< nl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        scalar L = 0;
        Info<< runTime.timeName();

        forAll(fieldNames, i)
        {
            const volScalarField base
            (
                IOobject
                (
                    fieldNames[i], runTime.timeName(), mesh,
                    IOobject::MUST_READ, IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh
            );
            const tmp<volScalarField> tpsi =
                conormalIterations > 0
              ? conormalPotential
                (
                    sub, cellMap, base, analyticXY, G, conormalIterations
                )
              : sealedPotential(sub, cellMap, base, analyticXY);

            tmp<volScalarField> tlap;
            if (insulated)
            {
                tlap = fvc::laplacian
                (
                    insulatedFaceConductivity(G, tpsi(), true),
                    tpsi()
                );
            }
            else
            {
                tlap = fvc::laplacian(G, tpsi());
            }

            if (args.found("insulationDelta"))
            {
                const tmp<volScalarField> tlapIns = fvc::laplacian
                (
                    insulatedFaceConductivity(G, tpsi(), true),
                    tpsi()
                );
                const tmp<volScalarField> tlapRaw = fvc::laplacian(G, tpsi());
                const scalarField d
                (
                    mag(tlapIns().primitiveField() - tlapRaw().primitiveField())
                );
                label iMax = findMax(d);
                Info<< " delta_" << fieldNames[i] << "=" << gMax(d)
                    << "/" << gMax(mag(tlapRaw().primitiveField()))
                    << "@cell" << iMax
                    << "(" << sub.C()[iMax].x() << ")";
            }

            const scalar Lpsi =
                gSum(tlap().primitiveField()*sub.V().field());
            const scalar LpsiNorm =
                gSum(mag(tlap().primitiveField())*sub.V().field());
            L += Lpsi;
            Info<< ' ' << Lpsi << ' ' << LpsiNorm;
        }
        Info<< ' ' << L;

        if (args.found("restDrift") && subsetterPtr && exposedI >= 0)
        {
            const volScalarField VmBase
            (
                IOobject
                (
                    "Vm", runTime.timeName(), mesh,
                    IOobject::MUST_READ, IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh
            );
            const scalar restValue(args.getOrDefault<scalar>("restValue", -0.084));
            boolList wall(sub.nCells(), false);
            for (const label cellI : sub.boundary()[exposedI].faceCells())
            {
                wall[cellI] = true;
            }
            scalar driftWall = 0, driftInner = 0;
            forAll(cellMap, cellI)
            {
                const scalar d = mag(VmBase[cellMap[cellI]] - restValue);
                if (wall[cellI])
                {
                    driftWall = max(driftWall, d);
                }
                else
                {
                    driftInner = max(driftInner, d);
                }
            }
            reduce(driftWall, maxOp<scalar>());
            reduce(driftInner, maxOp<scalar>());
            Info<< ' ' << driftWall << ' ' << driftInner;
        }

        if (args.found("gradientBias") && subsetterPtr && exposedI >= 0)
        {
            const volScalarField phiEBase
            (
                IOobject
                (
                    "phiE", runTime.timeName(), mesh,
                    IOobject::MUST_READ, IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh
            );
            const volVectorField gradGlobal(fvc::grad(phiEBase));
            const word trace
            (
                args.getOrDefault<word>("phiETrace", "zeroGradient")
            );
            autoPtr<volSymmTensorField> GePtr;
            autoPtr<volScalarField> sigmaBathPtr;
            if (trace == "conductivityWeighted")
            {
                GePtr.reset
                (
                    new volSymmTensorField
                    (
                        IOobject
                        (
                            args.getOrDefault<word>
                            (
                                "conductivityExtracellular",
                                "ConductivityExtracellular"
                            ),
                            "0", mesh, IOobject::MUST_READ,
                            IOobject::NO_WRITE, IOobject::NO_REGISTER
                        ),
                        mesh
                    )
                );
                sigmaBathPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            args.getOrDefault<word>
                            (
                                "sigmaBath", "bodyAndOrgansConductivity"
                            ),
                            "0", mesh, IOobject::MUST_READ,
                            IOobject::NO_WRITE, IOobject::NO_REGISTER
                        ),
                        mesh
                    )
                );
            }
            const tmp<volScalarField> tphiSub =
                trace == "zeroGradient"
              ? sealedPotential(sub, cellMap, phiEBase, false)
              : globalTracePotential
                (
                    sub, cellMap, subsetterPtr->faceMap(), phiEBase,
                    GePtr.get(), sigmaBathPtr.get()
                );
            const volVectorField gradSub(fvc::grad(tphiSub()));

            boolList wall(sub.nCells(), false);
            for (const label cellI : sub.boundary()[exposedI].faceCells())
            {
                wall[cellI] = true;
            }

            vectorField gRef(gradGlobal.primitiveField());
            const bool exactRef =
                args.getOrDefault<word>("exactPhiE", "") == "bathFDA";
            if (exactRef)
            {
                const scalar k(args.getOrDefault<scalar>("mmsK", 1/Foam::sqrt(2.0)));
                const scalar alpha(args.getOrDefault<scalar>("mmsAlpha", 0.01));
                const scalar se(args.get<scalar>("mmsSe"));
                const scalar s1t = Foam::sqrt(1 + runTime.value());
                const scalar pi = constant::mathematical::pi;
                forAll(gRef, cellI)
                {
                    const scalar x = mesh.C()[cellI].x();
                    gRef[cellI] =
                        vector
                        (
                            x <= 0 || x >= 1
                          ? 2*alpha/se
                          : k*s1t*pi*Foam::sin(pi*x) + alpha/se,
                            0,
                            0
                        );
                }
            }

            scalar wallDiff = 0, wallRef = 0, wallMax = 0;
            scalar innerDiff = 0, innerRef = 0, globalWallDiff = 0;
            forAll(gradSub, cellI)
            {
                const vector gG = gRef[cellMap[cellI]];
                const scalar d = mag(gradSub[cellI] - gG);
                if (wall[cellI])
                {
                    globalWallDiff +=
                        mag(gradGlobal[cellMap[cellI]] - gG);
                    wallDiff += d;
                    wallRef += mag(gG);
                    wallMax = max(wallMax, d/max(mag(gG), VSMALL));
                }
                else
                {
                    innerDiff += d;
                    innerRef += mag(gG);
                }
            }
            reduce(wallDiff, sumOp<scalar>());
            reduce(wallRef, sumOp<scalar>());
            reduce(wallMax, maxOp<scalar>());
            reduce(innerDiff, sumOp<scalar>());
            reduce(innerRef, sumOp<scalar>());
            reduce(globalWallDiff, sumOp<scalar>());
            Info<< ' ' << wallDiff/max(wallRef, VSMALL)
                << ' ' << wallMax
                << ' ' << innerDiff/max(innerRef, VSMALL);
            if (exactRef)
            {
                Info<< ' ' << globalWallDiff/max(wallRef, VSMALL);
            }
        }

        if (haveRef)
        {
            const volScalarField phiE
            (
                IOobject
                (
                    "phiE", runTime.timeName(), mesh,
                    IOobject::MUST_READ, IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh
            );
            if (refCell < 0)
            {
                refCell = mesh.findCell(args.get<point>("refPoint"));
                if (refCell < 0)
                {
                    FatalErrorInFunction
                        << "refPoint outside the mesh" << exit(FatalError);
                }
                const volScalarField sigmaBath
                (
                    IOobject
                    (
                        args.getOrDefault<word>
                        (
                            "sigmaBath", "bodyAndOrgansConductivity"
                        ),
                        "0", mesh, IOobject::MUST_READ, IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    mesh
                );
                const fvScalarMatrix Lap
                (
                    fvm::laplacian(fvc::interpolate(sigmaBath), phiE)
                );
                dRef = Lap.diag()[refCell];
            }
            const scalar predicted = refValue - L/dRef;
            const scalar solved = phiE[refCell];
            Info<< ' ' << dRef << ' ' << predicted << ' ' << solved << ' '
                << mag(predicted - solved)/max(mag(solved), VSMALL);
        }
        Info<< nl;
    }

    Info<< "End" << endl;
    return 0;
}


// ************************************************************************* //
