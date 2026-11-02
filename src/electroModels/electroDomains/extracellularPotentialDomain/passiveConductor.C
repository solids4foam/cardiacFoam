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

#include "passiveConductor.H"
#include "extracellularFaceConductivity.H"
#include "insulatedFaceConductivity.H"
#include "interfaceGradient.H"
#include "conormalStabilisation.H"
#include "PstreamReduceOps.H"
#include "fixedValueFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"

namespace Foam
{

const word passiveConductor::potentialName("phiE");


const Enum<passiveConductor::interfaceInterpolation>
passiveConductor::interfaceInterpolationNames
({
    {
        interfaceInterpolation::distanceWeightedHarmonic,
        "distanceWeightedHarmonic"
    },
    { interfaceInterpolation::conormalHarmonic, "conormalHarmonic" }
});


wordList passiveConductor::gradientNames()
{
    return wordList
    ({
        "grad(" + potentialName + ")",
        "grad(" + potentialName + ".component(0))"
    });
}


passiveConductor::torsoSettings passiveConductor::torsoSettings::read
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    torsoSettings s;

    s.interface = interfaceInterpolationNames.get
    (
        "interfaceConductivityInterpolation", dict
    );
    s.cellZones = dict.get<wordList>("cellZones");
    s.conductivityField = dict.get<word>("conductivityField");

    const dictionary* groundDict = dict.findDict("groundPatches");
    const dictionary* currentDict = dict.findDict("surfaceCurrentPatches");

    if (currentDict)
    {
        s.surfaceCurrentPatches = currentDict->toc();
        s.surfaceCurrentValues.setSize(s.surfaceCurrentPatches.size(), 0.0);

        forAll(s.surfaceCurrentPatches, patchI)
        {
            s.surfaceCurrentValues[patchI] =
                currentDict->get<scalar>(s.surfaceCurrentPatches[patchI]);
        }
    }

    if (groundDict)
    {
        s.groundPatches = groundDict->toc();

        forAll(s.groundPatches, patchI)
        {
            if (mesh.boundaryMesh().findPatchID(s.groundPatches[patchI]) < 0)
            {
                FatalErrorInFunction
                    << "Cannot find torso ground patch '"
                    << s.groundPatches[patchI] << "' on mesh '" << mesh.name()
                    << "'."
                    << exit(FatalError);
            }
        }
    }

    if (s.pinsReference())
    {
        if (!dict.found("phiERefPoint"))
        {
            FatalIOErrorInFunction(dict)
                << "The torso has no groundPatches, so phiERefPoint is "
                << "required: the potential is pinned at the cell that "
                << "contains it. Add phiERefPoint to the torso block, or a "
                << "ground patch."
                << exit(FatalIOError);
        }

        s.phiERefPoint = dict.get<point>("phiERefPoint");
        s.phiEReferenceValue = dict.get<scalar>("phiEReferenceValue");
    }
    else
    {
        // Ground patches fix the constant: no cell is pinned and the
        // reference is not read
        s.phiERefPoint = Zero;
        s.phiEReferenceValue = 0;
    }

    forAll(s.surfaceCurrentPatches, patchI)
    {
        if (mesh.boundaryMesh().findPatchID(s.surfaceCurrentPatches[patchI]) < 0)
        {
            FatalErrorInFunction
                << "Cannot find torso surface-current "
                << "patch '" << s.surfaceCurrentPatches[patchI]
                << "' on mesh '" << mesh.name() << "'."
                << exit(FatalError);
        }

        if (s.groundPatches.found(s.surfaceCurrentPatches[patchI]))
        {
            FatalErrorInFunction
                << "Torso patch '"
                << s.surfaceCurrentPatches[patchI]
                << "' cannot be listed in both groundPatches and "
                << "surfaceCurrentPatches."
                << exit(FatalError);
        }
    }

    if (groundDict)
    {
        s.groundValues.setSize(s.groundPatches.size(), 0.0);

        forAll(s.groundPatches, patchI)
        {
            s.groundValues[patchI] =
                groundDict->get<scalar>(s.groundPatches[patchI]);
        }
    }

    return s;
}


wordList passiveConductor::torsoSettings::patchTypes
(
    const fvMesh& mesh,
    const bool surfaceCurrents
) const
{
    wordList types
    (
        mesh.boundary().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (const word& patch : groundPatches)
    {
        types[mesh.boundaryMesh().findPatchID(patch)] =
            fixedValueFvPatchScalarField::typeName;
    }

    if (surfaceCurrents)
    {
        for (const word& patch : surfaceCurrentPatches)
        {
            types[mesh.boundaryMesh().findPatchID(patch)] =
                fixedGradientFvPatchScalarField::typeName;
        }
    }

    return types;
}


void passiveConductor::torsoSettings::setGroundValues
(
    volScalarField& psi
) const
{
    forAll(groundPatches, i)
    {
        const label patchI =
            psi.mesh().boundaryMesh().findPatchID(groundPatches[i]);

        psi.boundaryFieldRef()[patchI] == groundValues[i];
    }
}


label passiveConductor::locateReferenceCell
(
    const fvMesh& mesh,
    const point& p
)
{
    const label refCell = mesh.findCell(p);
    const label ownerCount =
        returnReduce(label(refCell >= 0 ? 1 : 0), sumOp<label>());

    if (ownerCount != 1)
    {
        FatalErrorInFunction
            << "phiERefPoint " << p << " must lie in exactly one cell of "
            << "one processor of mesh '" << mesh.name() << "', but "
            << ownerCount << " processors reported a containing cell."
            << exit(FatalError);
    }

    return refCell;
}


passiveConductor::passiveConductor
(
    const fvMesh& baseMesh,
    const fvMesh& heartMesh,
    const labelUList& heartCellMap,
    const labelUList& heartFaceMap,
    const volTensorField& Gi,
    const volTensorField& Ge,
    const torsoSettings& torso,
    const Switch sealedHeartBoundary,
    const volScalarField& sealingField,
    const bool totalWallConductivity,
    const bool registerFields
)
:
    baseMesh_(baseMesh),
    heartMesh_(heartMesh),
    Gi_(Gi),
    sealingField_(sealingField),
    torso_(torso),
    sealedHeartBoundary_(sealedHeartBoundary),
    nInterfaceGradientCorrectors_(0),
    totalWallConductivity_(totalWallConductivity),
    registerFields_(registerFields),
    heartCellToBaseCell_(),
    heartFaceToBaseFace_(),
    sigmaTotalPtr_(),
    sigmaIglobalPtr_(),
    sigmaExtracellularfPtr_(),
    referenceCell_(-1)
{
    if (torso_.interface == interfaceInterpolation::conormalHarmonic)
    {
        const dictionary& solverDict = baseMesh_.solverDict(potentialName);
        nInterfaceGradientCorrectors_ =
            solverDict.get<label>("nInterfaceGradientCorrectors");

        if (nInterfaceGradientCorrectors_ < 1)
        {
            FatalIOErrorInFunction(solverDict)
                << "nInterfaceGradientCorrectors must be at least 1; found "
                << nInterfaceGradientCorrectors_ << '.' << exit(FatalIOError);
        }
    }

    const dimensionSet conductivityDim
    (
        pow3(dimTime) * sqr(dimCurrent)/(dimMass*dimVolume)
    );

    sigmaTotalPtr_.reset
    (
        new volTensorField
        (
            IOobject
            (
                "sigmaTotal",
                baseMesh_.time().timeName(),
                baseMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                registerFields_
            ),
            baseMesh_,
            dimensionedTensor("zero", conductivityDim, tensor::zero),
            "zeroGradient"
        )
    );

    sigmaIglobalPtr_.reset
    (
        new volTensorField
        (
            IOobject
            (
                "sigmaIglobal",
                baseMesh_.time().timeName(),
                baseMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                registerFields_
            ),
            baseMesh_,
            dimensionedTensor("zero", conductivityDim, tensor::zero),
            "zeroGradient"
        )
    );

    buildHeartScatterMap(heartCellMap, heartFaceMap);

    // The torso conductivity field at its latest instance
    const volScalarField torsoConductivity
    (
        IOobject
        (
            torso_.conductivityField,
            baseMesh_.time().findInstance
            (
                baseMesh_.dbDir(), torso_.conductivityField
            ),
            baseMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            registerFields_
        ),
        baseMesh_
    );
    assembleConductivities(Ge, torsoConductivity);

    sigmaExtracellularfPtr_.reset
    (
        faceConductivity
        (
            sigmaTotalPtr_(),
            sigmaIglobalPtr_(),
            heartCellToBaseCell_,
            totalWallConductivity_,
            registerFields_
        ).ptr()
    );

    if (torso_.pinsReference())
    {
        referenceCell_ = locateReferenceCell(baseMesh_, torso_.phiERefPoint);
    }
}


void passiveConductor::buildHeartScatterMap
(
    const labelUList& heartCellMap,
    const labelUList& heartFaceMap
)
{
    const fvMesh& heartMesh = heartMesh_;

    heartCellToBaseCell_ = labelList(heartCellMap);

    if (heartCellToBaseCell_.size() != heartMesh.nCells())
    {
        FatalErrorInFunction
            << "Heart/base cell map size " << heartCellToBaseCell_.size()
            << " does not match myocardium cell count "
            << heartMesh.nCells() << "."
            << exit(FatalError);
    }

    labelList mappedHeartCell(baseMesh_.nCells(), -1);
    forAll(heartCellToBaseCell_, heartCellI)
    {
        const label baseCellI = heartCellToBaseCell_[heartCellI];

        if (baseCellI < 0 || baseCellI >= baseMesh_.nCells())
        {
            FatalErrorInFunction
                << "Heart cell " << heartCellI << " maps to invalid base cell "
                << baseCellI << "."
                << exit(FatalError);
        }

        if (mappedHeartCell[baseCellI] >= 0)
        {
            FatalErrorInFunction
                << "Heart cells " << mappedHeartCell[baseCellI] << " and "
                << heartCellI << " both map to base cell " << baseCellI << "."
                << exit(FatalError);
        }

        mappedHeartCell[baseCellI] = heartCellI;
    }

    {
        if (heartFaceMap.size() != heartMesh.nFaces())
        {
            FatalErrorInFunction
                << "Heart/base face map size " << heartFaceMap.size()
                << " does not match myocardium face count "
                << heartMesh.nFaces() << "."
                << exit(FatalError);
        }

        forAll(heartFaceMap, heartFaceI)
        {
            const label baseFaceI = heartFaceMap[heartFaceI];
            if (baseFaceI < 0 || baseFaceI >= baseMesh_.nFaces())
            {
                FatalErrorInFunction
                    << "Heart face " << heartFaceI
                    << " maps to invalid base face " << baseFaceI << "."
                    << exit(FatalError);
            }
        }

        heartFaceToBaseFace_ = labelList(heartFaceMap);
    }
}

void passiveConductor::assembleConductivities
(
    const volTensorField& Ge,
    const volScalarField& torsoConductivity
)
{
    volTensorField& sigmaTotal = sigmaTotalPtr_();
    volTensorField& sigmaIglobal = sigmaIglobalPtr_();
    tensorField& sigmaTotalI = sigmaTotal.primitiveFieldRef();
    tensorField& sigmaIglobalI = sigmaIglobal.primitiveFieldRef();
    const tensorField& GiI = Gi_.primitiveField();
    const tensorField& GeI = Ge.primitiveField();

    if
    (
        heartCellToBaseCell_.size() != GiI.size()
     || heartCellToBaseCell_.size() != GeI.size()
    )
    {
        FatalErrorInFunction
            << "Heart/base cell map size " << heartCellToBaseCell_.size()
            << " does not match intracellular/extracellular conductivity "
            << "cell counts " << GiI.size() << "/" << GeI.size() << "."
            << exit(FatalError);
    }

    labelList materialRegion(baseMesh_.nCells(), 0);

    forAll(heartCellToBaseCell_, heartCellI)
    {
        const label baseCellI = heartCellToBaseCell_[heartCellI];
        const tensor heartSigma = GiI[heartCellI] + GeI[heartCellI];

        if (magSqr(heartSigma) <= SMALL)
        {
            FatalErrorInFunction
                << "Myocardium cell " << heartCellI
                << " has zero total conductivity."
                << exit(FatalError);
        }

        sigmaTotalI[baseCellI] = heartSigma;
        sigmaIglobalI[baseCellI] = GiI[heartCellI];
        materialRegion[baseCellI] = 1;
    }

    const scalarField& bathSigmaI = torsoConductivity.primitiveField();

    forAll(torso_.cellZones, zoneNameI)
    {
        const word& zoneName = torso_.cellZones[zoneNameI];
        const label zoneI = baseMesh_.cellZones().findZoneID(zoneName);

        if (zoneI < 0)
        {
            FatalErrorInFunction
                << "Bath cellZone '" << zoneName
                << "' not found on base mesh '" << baseMesh_.name() << "'."
                << exit(FatalError);
        }

        const labelList& cells = baseMesh_.cellZones()[zoneI];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const scalar sigma = bathSigmaI[cellI];

            if (materialRegion[cellI] != 0)
            {
                FatalErrorInFunction
                    << "Base cell " << cellI << " in bath cellZone '"
                    << zoneName << "' is already assigned to "
                    <<
                    (
                        materialRegion[cellI] == 1
                      ? "myocardium"
                      : "another bath zone"
                    )
                    << "."
                    << exit(FatalError);
            }

            if (sigma <= SMALL)
            {
                FatalErrorInFunction
                    << "Bath cell " << cellI << " in cellZone '" << zoneName
                    << "' has non-positive conductivity " << sigma << "."
                    << exit(FatalError);
            }

            sigmaTotalI[cellI] = sigma*tensor::I;
            materialRegion[cellI] = 2;
        }
    }

    label unassignedCellCount = 0;
    forAll(materialRegion, cellI)
    {
        if (materialRegion[cellI] == 0)
        {
            ++unassignedCellCount;
        }
    }
    reduce(unassignedCellCount, sumOp<label>());

    if (unassignedCellCount != 0)
    {
        FatalErrorInFunction
            << unassignedCellCount << " base-mesh cells are assigned to "
            << "neither the myocardium nor a configured bath cellZone."
            << exit(FatalError);
    }

    sigmaTotal.correctBoundaryConditions();
    sigmaIglobal.correctBoundaryConditions();
}

namespace
{

// Face value between two torso cells: the distance-weighted harmonic mean
// of the isotropic zone conductivities, linear where they are equal.
tensor torsoFaceConductivity
(
    const tensor& a,
    const tensor& b,
    const scalar weightA
)
{
    if (a.xx() == b.xx())
    {
        return extracellularFaceConductivity::linear(a, b, weightA);
    }

    return
        extracellularFaceConductivity::distanceWeightedHarmonic
        (
            a.xx(), b.xx(), 1.0 - weightA, weightA
        )*tensor::I;
}


// Face conductivity between the cells a and b of a face with unit normal n
// and weight weightA of a: from their total and intracellular
// conductivities and whether each is a heart cell. A heart cell enters
// with Ge, or with Gi + Ge on a heart-torso face when totalWallConductivity.
tensor faceValue
(
    const tensor& sigmaA,
    const tensor& sigmaIA,
    const bool heartA,
    const tensor& sigmaB,
    const tensor& sigmaIB,
    const bool heartB,
    const vector& n,
    const scalar weightA,
    const bool totalWallConductivity
)
{
    const bool wallFace = heartA != heartB;
    const tensor sigmaEA =
        heartA && !(wallFace && totalWallConductivity)
      ? sigmaA - sigmaIA
      : sigmaA;
    const tensor sigmaEB =
        heartB && !(wallFace && totalWallConductivity)
      ? sigmaB - sigmaIB
      : sigmaB;

    if (!wallFace)
    {
        return
            heartA
          ? extracellularFaceConductivity::linear(sigmaEA, sigmaEB, weightA)
          : torsoFaceConductivity(sigmaEA, sigmaEB, weightA);
    }

    return extracellularFaceConductivity::conormalHarmonic
    (
        sigmaEA, sigmaEB, n, 1.0 - weightA, weightA
    );
}

} // End anonymous namespace


tmp<surfaceTensorField> passiveConductor::faceConductivity
(
    const volTensorField& sigma,
    const volTensorField& sigmaIglobal,
    const labelUList& heartCells,
    const bool totalWallConductivity,
    const bool registerObject
)
{
    const fvMesh& m = sigma.mesh();

    tmp<surfaceTensorField> tSef
    (
        new surfaceTensorField
        (
            IOobject
            (
                "sigmaExtracellularf",
                m.time().timeName(),
                m,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                registerObject
            ),
            m,
            dimensionedTensor("zero", sigma.dimensions(), tensor::zero)
        )
    );

    surfaceTensorField& Sef = tSef.ref();

    // Myocardium membership from the heart cell map, with coupled neighbours
    volScalarField heartFlag
    (
        IOobject
        (
            "heartFlag",
            m.time().timeName(),
            m,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        m,
        dimensionedScalar(dimless, 0)
    );
    for (const label baseCellI : heartCells)
    {
        heartFlag[baseCellI] = 1;
    }
    heartFlag.correctBoundaryConditions();

    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const tensorField& sigmaI = sigma.primitiveField();
    const tensorField& sigmaIntracellularI = sigmaIglobal.primitiveField();
    const scalarField& weights = m.weights().primitiveField();
    tensorField& SefI = Sef.primitiveFieldRef();

    forAll(own, faceI)
    {
        const label ownCell = own[faceI];
        const label neiCell = nei[faceI];

        SefI[faceI] = faceValue
        (
            sigmaI[ownCell],
            sigmaIntracellularI[ownCell],
            heartFlag[ownCell] > 0.5,
            sigmaI[neiCell],
            sigmaIntracellularI[neiCell],
            heartFlag[neiCell] > 0.5,
            m.Sf()[faceI]/m.magSf()[faceI],
            weights[faceI],
            totalWallConductivity
        );
    }

    forAll(Sef.boundaryField(), patchI)
    {
        if (sigma.boundaryField()[patchI].coupled())
        {
            const tensorField sigmaP
            (
                sigma.boundaryField()[patchI].patchInternalField()
            );
            const tensorField sigmaN
            (
                sigma.boundaryField()[patchI].patchNeighbourField()
            );
            const tensorField sigmaIP
            (
                sigmaIglobal.boundaryField()[patchI].patchInternalField()
            );
            const tensorField sigmaIN
            (
                sigmaIglobal.boundaryField()[patchI].patchNeighbourField()
            );
            const scalarField heartP
            (
                heartFlag.boundaryField()[patchI].patchInternalField()
            );
            const scalarField heartN
            (
                heartFlag.boundaryField()[patchI].patchNeighbourField()
            );

            Field<tensor>& Sefp = Sef.boundaryFieldRef()[patchI];
            const scalarField& patchWeights =
                m.weights().boundaryField()[patchI];
            const vectorField patchNf(m.boundary()[patchI].nf());

            forAll(Sefp, faceI)
            {
                Sefp[faceI] = faceValue
                (
                    sigmaP[faceI],
                    sigmaIP[faceI],
                    heartP[faceI] > 0.5,
                    sigmaN[faceI],
                    sigmaIN[faceI],
                    heartN[faceI] > 0.5,
                    patchNf[faceI],
                    patchWeights[faceI],
                    totalWallConductivity
                );
            }
        }
        else
        {
            // Physical boundaries use the cell-internal extracellular value;
            // the heart block carries Gi there.
            Sef.boundaryFieldRef()[patchI] =
                sigma.boundaryField()[patchI].patchInternalField()
              - sigmaIglobal.boundaryField()[patchI].patchInternalField();
        }
    }

    return tSef;
}


tmp<volVectorField> passiveConductor::jumpAwareGradient
(
    const volScalarField& psi
) const
{
    // Heart-side wall conductivity for the two-sided flux continuity
    const volTensorField sigmaWall
    (
        totalWallConductivity_
      ? volTensorField(sigmaTotalPtr_())
      : volTensorField(sigmaTotalPtr_() - sigmaIglobalPtr_())
    );
    tmp<volVectorField> tgrad(fvc::grad(psi));
    for (label iter = 0; iter < nInterfaceGradientCorrectors_; ++iter)
    {
        tgrad = interfaceGradient
        (
            psi, sigmaWall, sigmaIglobalPtr_(), tgrad()
        );
    }
    return tgrad;
}


tmp<fvScalarMatrix> passiveConductor::conormalOperator
(
    volScalarField& psi,
    const volScalarField* source,
    faceCorrection* correction
) const
{
    const surfaceTensorField& Sef = sigmaExtracellularfPtr_();
    const surfaceVectorField& Sf = baseMesh_.Sf();
    const surfaceScalarField gammaNN
    (
        (Sf & Sef & Sf)/sqr(baseMesh_.magSf())
    );
    const surfaceVectorField SfGammaCorr((Sf & Sef) - gammaNN*Sf);
    tmp<volVectorField> tgrad(jumpAwareGradient(psi));

    tmp<fvScalarMatrix> tEqn
    (
        fvm::laplacian(gammaNN, psi)
      + conormalStabilisation(Sef, psi)
      + fvc::div(SfGammaCorr & linearInterpolate(tgrad()))
    );
    if (source)
    {
        tEqn.ref() -= *source;
    }

    // Non-orthogonal correction of fvm::laplacian rebuilt from the
    // jump-aware gradient
    const surfaceScalarField nonOrthGradCorr
    (
        gammaNN*baseMesh_.magSf()
       *(
            baseMesh_.nonOrthCorrectionVectors()
          & linearInterpolate(tgrad() - fvc::grad(psi))
        )
    );
    tEqn.ref() += fvc::div(nonOrthGradCorr);

    if (correction)
    {
        correction->crossFlux = SfGammaCorr & linearInterpolate(tgrad());
        correction->crossFlux.ref() += nonOrthGradCorr;
        correction->gradient = tgrad;
    }

    return tEqn;
}


tmp<fvScalarMatrix> passiveConductor::harmonicOperator
(
    volScalarField& psi,
    const volScalarField* source,
    faceCorrection* correction
) const
{
    tmp<fvScalarMatrix> tEqn
    (
        fvm::laplacian(sigmaExtracellularfPtr_(), psi)
      + conormalStabilisation(sigmaExtracellularfPtr_(), psi)
    );
    if (source)
    {
        tEqn.ref() -= *source;
    }

    if (correction)
    {
        correction->gradient = fvc::grad(psi);
    }

    return tEqn;
}


void passiveConductor::addHeartBlock
(
    fvScalarMatrix& eqn,
    volScalarField& heartPsi
) const
{
    const tmp<surfaceTensorField> tGif
    (
        insulatedFaceConductivity(Gi_, sealingField_, sealedHeartBoundary_)
    );
    fvScalarMatrix heartEqn
    (
        fvm::laplacian(tGif(), heartPsi)
      + conormalStabilisation(tGif(), heartPsi)
    );
    const labelUList& heartFaceMap = heartFaceToBaseFace_;
    const labelUList& heartOwner = heartPsi.mesh().owner();

    forAll(heartCellToBaseCell_, heartCellI)
    {
        const label baseCellI = heartCellToBaseCell_[heartCellI];
        eqn.diag()[baseCellI] += heartEqn.diag()[heartCellI];
        eqn.source()[baseCellI] += heartEqn.source()[heartCellI];
    }
    forAll(heartOwner, heartFaceI)
    {
        const label baseFaceI = heartFaceMap[heartFaceI];
        eqn.upper()[baseFaceI] += heartEqn.upper()[heartFaceI];
    }

    forAll(heartEqn.internalCoeffs(), patchI)
    {
        const fvPatch& heartPatch = heartPsi.mesh().boundary()[patchI];
        forAll(heartEqn.internalCoeffs()[patchI], faceI)
        {
            const scalar internalCoeff =
                heartEqn.internalCoeffs()[patchI][faceI];
            const scalar boundaryCoeff =
                heartEqn.boundaryCoeffs()[patchI][faceI];
            if
            (
                mag(internalCoeff) <= SMALL
             && mag(boundaryCoeff) <= SMALL
            )
            {
                continue;
            }

            const label heartFaceI = heartPatch.start() + faceI;
            const label baseFaceI = heartFaceMap[heartFaceI];
            if (baseFaceI < baseMesh_.nInternalFaces())
            {
                FatalErrorInFunction
                    << "Non-zero heart boundary coefficient maps to base "
                    << "internal face " << baseFaceI << ". Only coupled or "
                    << "physical base-boundary coefficients can be mapped."
                    << exit(FatalError);
            }

            const label basePatchI =
                baseMesh_.boundaryMesh().whichPatch(baseFaceI);
            const label basePatchFaceI =
                baseFaceI - baseMesh_.boundary()[basePatchI].start();
            eqn.internalCoeffs()[basePatchI][basePatchFaceI] +=
                internalCoeff;
            eqn.boundaryCoeffs()[basePatchI][basePatchFaceI] +=
                boundaryCoeff;
        }
    }
}


tmp<fvScalarMatrix> passiveConductor::assemble
(
    volScalarField& psi,
    volScalarField& heartPsi,
    const volScalarField* source,
    faceCorrection* correction
) const
{
    if (psi.name() != potentialName || heartPsi.name() != potentialName)
    {
        FatalErrorInFunction
            << "The conductor assembles the potential " << potentialName
            << ", not " << psi.name() << " and " << heartPsi.name() << '.'
            << exit(FatalError);
    }

    tmp<fvScalarMatrix> tEqn
    (
        torso_.interface == interfaceInterpolation::conormalHarmonic
      ? conormalOperator(psi, source, correction)
      : harmonicOperator(psi, source, correction)
    );
    addHeartBlock(tEqn.ref(), heartPsi);

    return tEqn;
}


void passiveConductor::fillHeart
(
    const volScalarField& psi,
    volScalarField& heartPsi
) const
{
    scalarField& heartPsiI = heartPsi.primitiveFieldRef();
    forAll(heartCellToBaseCell_, heartCellI)
    {
        heartPsiI[heartCellI] = psi[heartCellToBaseCell_[heartCellI]];
    }
    heartPsi.correctBoundaryConditions();
}


void passiveConductor::setSurfaceCurrents(volScalarField& psi) const
{
    const volTensorField& sigma = sigmaTotalPtr_();

    forAll(torso_.surfaceCurrentPatches, i)
    {
        const label patchI =
            baseMesh_.boundaryMesh().findPatchID
            (
                torso_.surfaceCurrentPatches[i]
            );

        fixedGradientFvPatchScalarField& psiP =
            refCast<fixedGradientFvPatchScalarField>
            (
                psi.boundaryFieldRef()[patchI]
            );

        const tensorField sigmaP
        (
            sigma.boundaryField()[patchI].patchInternalField()
        );
        const vectorField nf(psiP.patch().nf());
        scalarField& gradient = psiP.gradient();

        forAll(gradient, faceI)
        {
            const scalar sigmaN = nf[faceI] & (sigmaP[faceI] & nf[faceI]);

            if (sigmaN <= SMALL)
            {
                FatalErrorInFunction
                    << "Surface-current patch '"
                    << torso_.surfaceCurrentPatches[i]
                    << "' found non-positive normal conductivity " << sigmaN
                    << " on face " << faceI << "."
                    << exit(FatalError);
            }

            gradient[faceI] = torso_.surfaceCurrentValues[i]/sigmaN;
        }

        psiP.evaluate();
    }
}


void passiveConductor::applyReference
(
    fvScalarMatrix& eqn,
    const scalar value
) const
{
    // The cell is -1 when ground patches fix the constant: a pin would then
    // over-determine the problem
    if (referenceCell_ >= 0)
    {
        eqn.setReference(referenceCell_, value, true);
    }
}


namespace
{

// True for a gradient whose stencil is the cell and its face neighbours
bool faceStencilGradient(const fvMesh& mesh, const word& name)
{
    ITstream is(mesh.gradScheme(name));
    const word scheme(is);
    if (scheme == "leastSquares")
    {
        return true;
    }
    if (scheme == "Gauss" && !is.eof())
    {
        const word interpolation(is);
        return interpolation == "linear";
    }
    return false;
}

} // End anonymous namespace


List<probedOperator::hop> passiveConductor::explicitStencil() const
{
    bool faceStencil = true;
    for (const fvMesh* meshPtr : {&baseMesh_, &heartMesh_})
    {
        for (const word& name : gradientNames())
        {
            faceStencil = faceStencil && faceStencilGradient(*meshPtr, name);
        }
    }

    // Divergence and face interpolation, jump-aware gradient passes, and
    // the gradient stencil
    List<probedOperator::hop> hops(1, probedOperator::hop::face);
    for (label passI = 0; passI < nInterfaceGradientCorrectors_; ++passI)
    {
        hops.append(probedOperator::hop::face);
    }
    hops.append
    (
        faceStencil ? probedOperator::hop::face : probedOperator::hop::point
    );
    return hops;
}


tmp<volScalarField> passiveConductor::heartSource
(
    const volScalarField& Vm
) const
{
    tmp<volScalarField> tHeartRhs =
        -fvc::laplacian
        (
            insulatedFaceConductivity(Gi_, Vm, sealedHeartBoundary_),
            Vm
        );
    const volScalarField& heartRhs = tHeartRhs();

    tmp<volScalarField> trhsGlobal
    (
        new volScalarField
        (
            IOobject
            (
                "phiESource",
                baseMesh_.time().timeName(),
                baseMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                registerFields_
            ),
            baseMesh_,
            dimensionedScalar("zero", heartRhs.dimensions(), 0.0),
            "zeroGradient"
        )
    );
    volScalarField& rhsGlobal = trhsGlobal.ref();

    scalarField& rhsGlobalI = rhsGlobal.primitiveFieldRef();
    const scalarField& heartRhsI = heartRhs.primitiveField();
    rhsGlobalI = 0.0;

    forAll(heartCellToBaseCell_, heartCellI)
    {
        rhsGlobalI[heartCellToBaseCell_[heartCellI]] = heartRhsI[heartCellI];
    }

    rhsGlobal.correctBoundaryConditions();

    return trhsGlobal;
}

} // End namespace Foam

// ************************************************************************* //
