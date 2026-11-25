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
#include "referenceCell.H"
#include "extracellularFaceConductivity.H"
#include "insulatedFaceConductivity.H"
#include "sealedWallPolicy.H"
#include "interfaceGradient.H"
#include "centredCPCCellToCellStencilObject.H"
#include "SVD.H"
#include "syncTools.H"
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

    // The current density of each patch as a function of time: a plain
    // number is a constant
    s.surfaceCurrents.setSize(s.surfaceCurrentPatches.size());
    forAll(s.surfaceCurrentPatches, patchI)
    {
        s.surfaceCurrents.set
        (
            patchI,
            Function1<scalar>::New
            (
                s.surfaceCurrentPatches[patchI], *currentDict, &mesh
            )
        );
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


bool passiveConductor::torsoSettings::constantSurfaceCurrents() const
{
    for (const Function1<scalar>& current : surfaceCurrents)
    {
        if (!current.constant())
        {
            return false;
        }
    }

    return true;
}


scalarList passiveConductor::torsoSettings::surfaceCurrentValues
(
    const scalar t
) const
{
    scalarList values(surfaceCurrents.size());

    forAll(surfaceCurrents, i)
    {
        values[i] = surfaceCurrents[i].value(t);
    }

    return values;
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


bool passiveConductor::totalWallConductivity
(
    const torsoCoupling::type coupling
)
{
    switch (coupling)
    {
        case torsoCoupling::type::bath:
            return false;
        case torsoCoupling::type::leadField:
            return true;
        case torsoCoupling::type::none:
            break;
    }

    FatalErrorInFunction
        << "A passive conductor serves the bath or the lead fields, not "
        << torsoCoupling::name(coupling) << '.' << exit(FatalError);

    return false;
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
    const torsoCoupling::type coupling,
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
    coupling_(coupling),
    totalWallConductivity_(totalWallConductivity(coupling)),
    registerFields_(registerFields),
    heartCellToBaseCell_(),
    heartFaceToBaseFace_(),
    sigmaTotalPtr_(),
    sigmaIglobalPtr_(),
    sigmaExtracellularfPtr_(),
    referenceCell_(-1),
    wallTransferBuilt_(false),
    wallCell_(),
    wallValueStart_(),
    wallValueIndex_(),
    wallGradientStart_(),
    wallGradientIndex_(),
    wallRowStart_(),
    internalLayerHeart_(),
    internalLayerTorso_(),
    internalLayerWall_(),
    coupledLayerHeart_(),
    coupledLayerSlot_(),
    coupledLayerWall_(),
    coupledReceiveCell_(),
    coupledReceiveSlot_(),
    layerCoeffStart_(),
    layerCoeff_()
{
    // Lead fields are built once and hold the surface currents only in the
    // constant term of each channel; the heart of a leadField coupling does
    // not see them. A current that varies in time has no representation.
    if
    (
        coupling_ == torsoCoupling::type::leadField
     && !torso_.constantSurfaceCurrents()
    )
    {
        wordList varying;
        forAll(torso_.surfaceCurrents, i)
        {
            if (!torso_.surfaceCurrents[i].constant())
            {
                varying.append(torso_.surfaceCurrentPatches[i]);
            }
        }

        FatalErrorInFunction
            << "The surface currents of patches " << flatOutput(varying)
            << " vary in time, but torsoCoupling leadField needs constant "
            << "ones: the lead fields are built once and carry the surface "
            << "currents only in the constant term of each ECG channel, and "
            << "the heart does not see them. Use torsoCoupling bath for a "
            << "current that varies in time."
            << exit(FatalError);
    }

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


void passiveConductor::setSurfaceCurrents
(
    volScalarField& psi,
    const scalar t
) const
{
    const volTensorField& sigma = sigmaTotalPtr_();

    forAll(torso_.surfaceCurrentPatches, i)
    {
        const scalar current = torso_.surfaceCurrents[i].value(t);

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

            gradient[faceI] = current/sigmaN;
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


// Mark the cells of mesh whose gradient name reads beyond the face
// neighbours: none for a face-stencil scheme, the cells with a face on a
// sealed or value patch for sealedWallLeastSquares (whose other cells take
// the face-neighbour least squares), every cell otherwise
void markPointStencilCells
(
    const fvMesh& mesh,
    const word& name,
    boolList& mark
)
{
    if (faceStencilGradient(mesh, name))
    {
        return;
    }

    ITstream is(mesh.gradScheme(name));
    const word scheme(is);
    if (scheme != "sealedWallLeastSquares")
    {
        mark = true;
        return;
    }

    const labelList& states =
        sealedWallPolicy::states(mesh, sealedWallPolicy::fieldKey(name));
    forAll(states, patchI)
    {
        if (states[patchI] != sealedWallPolicy::stock)
        {
            for (const label cellI : mesh.boundary()[patchI].faceCells())
            {
                mark[cellI] = true;
            }
        }
    }
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


boolList passiveConductor::explicitStencilPointCells() const
{
    boolList base(baseMesh_.nCells(), false);
    boolList heart(heartMesh_.nCells(), false);
    for (const word& name : gradientNames())
    {
        markPointStencilCells(baseMesh_, name, base);
        markPointStencilCells(heartMesh_, name, heart);
    }
    forAll(heart, heartCellI)
    {
        if (heart[heartCellI])
        {
            base[heartCellToBaseCell_[heartCellI]] = true;
        }
    }
    return base;
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

    if (totalWallConductivity_)
    {
        transferWallLayer(Vm, rhsGlobal);
        rhsGlobal.correctBoundaryConditions();
    }

    return trhsGlobal;
}


namespace
{

// 1 in the heart cells of the base mesh, 0 in the torso; coupled patches
// carry the neighbour values
tmp<volScalarField> heartIndicator
(
    const fvMesh& m,
    const labelUList& heartCells
)
{
    tmp<volScalarField> tFlag
    (
        new volScalarField
        (
            IOobject
            (
                "heartIndicator",
                m.time().timeName(),
                m,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            m,
            dimensionedScalar(dimless, 0)
        )
    );
    for (const label baseCellI : heartCells)
    {
        tFlag.ref()[baseCellI] = 1;
    }
    tFlag.ref().correctBoundaryConditions();

    return tFlag;
}


// Torso share w = kT/(kH + kT), k = n.sigma.n/d, of a wall layer on the
// face between a heart and a torso cell at distance shares dH and dT
scalar torsoShare
(
    const tensor& sigmaHeart,
    const tensor& sigmaTorso,
    const vector& n,
    const scalar dH,
    const scalar dT
)
{
    const scalar kH = (n & sigmaHeart & n)/max(dH, SMALL);
    const scalar kT = (n & sigmaTorso & n)/max(dT, SMALL);

    return kT/max(kH + kT, VSMALL);
}

} // End anonymous namespace


namespace
{

// Least-squares quadratic of Vm about a heart wall cell: the stencil
// entries of its value and gradient rows, the scale 1/(L|d|) of each value
// row, the length L and the pseudo-inverse of the row matrix
struct wallFit
{
    scalar L;
    labelList valueEntries;
    scalarList valueScale;
    labelList gradientEntries;
    scalarRectangularMatrix pinv;
};


// fld in the flat (compact) addressing of map: internal values, boundary
// values, then the remote values the stencil reads; as the first step of
// extendedCellToFaceStencil::collectData
template<class Type>
List<Type> flatStencilField
(
    const mapDistribute& map,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    List<Type> flat(map.constructSize(), Zero);

    forAll(fld, celli)
    {
        flat[celli] = fld[celli];
    }
    forAll(fld.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchi];
        label nCompact =
            pfld.patch().start()
          - fld.mesh().nInternalFaces()
          + fld.mesh().nCells();

        forAll(pfld, i)
        {
            flat[nCompact++] = pfld[i];
        }
    }

    map.distribute(flat);

    return flat;
}

} // End anonymous namespace


void passiveConductor::buildWallTransfer(const fvMesh& heartMesh) const
{
    if (heartMesh.changing())
    {
        FatalErrorInFunction
            << "The wall layer is built once; the heart mesh "
            << heartMesh.name() << " changes." << exit(FatalError);
    }

    const vectorField& heartC = heartMesh.C().primitiveField();

    // Whether a stencil entry is a cell (not a boundary face), and whether
    // it is a cell without non-coupled boundary faces, whose gradient is
    // not one-sided
    const extendedCentredCellToCellStencil& stencil =
        centredCPCCellToCellStencilObject::New(heartMesh);
    volScalarField isCell
    (
        IOobject
        (
            "wallLayerIsCell",
            heartMesh.time().timeName(),
            heartMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        heartMesh,
        dimensionedScalar(dimless, 1)
    );
    volScalarField isInterior("wallLayerIsInterior", isCell);
    forAll(heartMesh.boundary(), patchI)
    {
        isCell.boundaryFieldRef()[patchI] = 0;
        isInterior.boundaryFieldRef()[patchI] = 0;
        if (!heartMesh.boundary()[patchI].coupled())
        {
            for (const label c : heartMesh.boundary()[patchI].faceCells())
            {
                isInterior[c] = 0;
            }
        }
    }
    List<List<vector>> stencilC;
    List<List<scalar>> stencilIsCell;
    List<List<scalar>> stencilIsInterior;
    stencil.collectData(heartMesh.C(), stencilC);
    stencil.collectData(isCell, stencilIsCell);
    stencil.collectData(isInterior, stencilIsInterior);

    // The heart wall cells: heart cells of the heart-torso faces
    const tmp<volScalarField> tHeart
    (
        heartIndicator(baseMesh_, heartCellToBaseCell_)
    );
    const volScalarField& heart = tHeart();
    labelList baseToHeart(baseMesh_.nCells(), -1);
    forAll(heartCellToBaseCell_, heartCellI)
    {
        baseToHeart[heartCellToBaseCell_[heartCellI]] = heartCellI;
    }
    boolList wall(heartMesh.nCells(), false);
    const labelUList& own = baseMesh_.owner();
    const labelUList& nei = baseMesh_.neighbour();
    forAll(own, faceI)
    {
        const bool ownHeart = heart[own[faceI]] > 0.5;
        if (ownHeart != (heart[nei[faceI]] > 0.5))
        {
            wall[baseToHeart[ownHeart ? own[faceI] : nei[faceI]]] = true;
        }
    }
    forAll(baseMesh_.boundary(), patchI)
    {
        const fvPatch& p = baseMesh_.boundary()[patchI];
        if (!p.coupled())
        {
            continue;
        }
        const scalarField heartN
        (
            heart.boundaryField()[patchI].patchNeighbourField()
        );
        forAll(p.faceCells(), faceI)
        {
            const label c = p.faceCells()[faceI];
            if (heart[c] > 0.5 && heartN[faceI] < 0.5)
            {
                wall[baseToHeart[c]] = true;
            }
        }
    }

    // Quadratic in every geometric direction: a and the symmetric H
    const label nD = heartMesh.nGeometricD();
    const label nUnknowns = nD + nD*(nD + 1)/2;

    labelList fitIndex(heartMesh.nCells(), -1);
    List<wallFit> fits(findIndices(wall, true).size());
    label nFits = 0;

    forAll(wall, c)
    {
        if (!wall[c])
        {
            continue;
        }
        const List<vector>& pts = stencilC[c];
        wallFit& fit = fits[nFits];
        fitIndex[c] = nFits++;

        fit.L = 0;
        forAll(pts, k)
        {
            if (stencilIsCell[c][k] > 0.5)
            {
                fit.L = max(fit.L, mag(pts[k] - heartC[c]));
            }
        }

        // Value rows (V_k - V_c)/(L|d|) = (a.d + d.H.d/2)/|d| and gradient
        // rows grad V_k = a + H d, with d scaled by L and H by L
        DynamicList<FixedList<scalar, 9>> rows(4*pts.size());
        DynamicList<label> valueEntries(pts.size());
        DynamicList<scalar> valueScale(pts.size());
        DynamicList<label> gradientEntries(pts.size());
        forAll(pts, k)
        {
            if (stencilIsCell[c][k] < 0.5)
            {
                continue;
            }
            const vector d((pts[k] - heartC[c])/fit.L);
            const scalar md = mag(d);
            if (md > SMALL)
            {
                FixedList<scalar, 9> r;
                r[0] = d.x()/md;
                r[1] = d.y()/md;
                r[2] = d.z()/md;
                r[3] = 0.5*d.x()*d.x()/md;
                r[4] = d.x()*d.y()/md;
                r[5] = d.x()*d.z()/md;
                r[6] = 0.5*d.y()*d.y()/md;
                r[7] = d.y()*d.z()/md;
                r[8] = 0.5*d.z()*d.z()/md;
                rows.append(r);
                valueEntries.append(k);
                valueScale.append(1/(fit.L*md));
            }
        }
        forAll(pts, k)
        {
            if (stencilIsCell[c][k] < 0.5 || stencilIsInterior[c][k] < 0.5)
            {
                continue;
            }
            const vector d((pts[k] - heartC[c])/fit.L);
            FixedList<scalar, 9> rx(Zero);
            FixedList<scalar, 9> ry(Zero);
            FixedList<scalar, 9> rz(Zero);
            rx[0] = 1; rx[3] = d.x(); rx[4] = d.y(); rx[5] = d.z();
            ry[1] = 1; ry[4] = d.x(); ry[6] = d.y(); ry[7] = d.z();
            rz[2] = 1; rz[5] = d.x(); rz[7] = d.y(); rz[8] = d.z();
            rows.append(rx);
            rows.append(ry);
            rows.append(rz);
            gradientEntries.append(k);
        }

        scalarRectangularMatrix A(rows.size(), 9);
        forAll(rows, i)
        {
            for (label j = 0; j < 9; ++j)
            {
                A(i, j) = rows[i][j];
            }
        }
        const SVD svd(A, 1e-10);
        if (9 - svd.nZeros() < nUnknowns)
        {
            FatalErrorInFunction
                << "The quadratic fit of the wall layer at heart cell " << c
                << " (centre " << heartC[c] << ") has rank "
                << 9 - svd.nZeros() << " of " << nUnknowns << ": its "
                << "cell-point-cell stencil (" << valueEntries.size()
                << " cells, " << gradientEntries.size() << " without "
                << "boundary faces) does not determine grad(Vm) at the "
                << "heart-torso wall." << exit(FatalError);
        }
        fit.pinv = svd.VSinvUt();
        fit.valueEntries.transfer(valueEntries);
        fit.valueScale.transfer(valueScale);
        fit.gradientEntries.transfer(gradientEntries);
    }

    // Row addressing of the wall cells, in flat stencil indices
    const labelListList& compact = stencil.stencil();
    wallCell_.setSize(nFits);
    wallValueStart_.setSize(nFits + 1);
    wallGradientStart_.setSize(nFits + 1);
    wallRowStart_.setSize(nFits + 1);
    label nValue = 0, nGradient = 0;
    forAll(fits, w)
    {
        nValue += fits[w].valueEntries.size();
        nGradient += fits[w].gradientEntries.size();
    }
    wallValueIndex_.setSize(nValue);
    wallGradientIndex_.setSize(nGradient);
    wallValueStart_[0] = 0;
    wallGradientStart_[0] = 0;
    wallRowStart_[0] = 0;
    forAll(wall, c)
    {
        const label w = fitIndex[c];
        if (w < 0)
        {
            continue;
        }
        const wallFit& fit = fits[w];
        wallCell_[w] = c;
        label iv = wallValueStart_[w];
        for (const label k : fit.valueEntries)
        {
            wallValueIndex_[iv++] = compact[c][k];
        }
        wallValueStart_[w + 1] = iv;
        label ig = wallGradientStart_[w];
        for (const label k : fit.gradientEntries)
        {
            wallGradientIndex_[ig++] = compact[c][k];
        }
        wallGradientStart_[w + 1] = ig;
        wallRowStart_[w + 1] =
            wallRowStart_[w]
          + fit.valueEntries.size() + 3*fit.gradientEntries.size();
    }

    // The layer of a face, w q|S| with q|S| = S.Gi.grad(Vm)(x) and
    // grad(Vm)(x) = J(e) u, u = pinv r, e = (x - x_c)/L, is linear in the
    // rows r of its wall cell: layer = (pinv^T J(e)^T a).r with
    // a = w Gi^T S. The value-row scales go into the coefficients.
    const volTensorField& sigma = sigmaTotalPtr_();
    const surfaceVectorField& Sf = baseMesh_.Sf();
    const surfaceVectorField& Cf = baseMesh_.Cf();
    const surfaceScalarField& magSf = baseMesh_.magSf();
    const scalarField& weights = baseMesh_.weights().primitiveField();

    DynamicList<scalar> coeff;
    DynamicList<label> coeffStart;
    coeffStart.append(0);

    const auto addLayer = [&]
    (
        const label h,
        const tensor& sigmaTorso,
        const vector& S,
        const scalar magS,
        const point& x,
        const scalar dH
    ) -> label
    {
        const label hHeart = baseToHeart[h];
        const label w = fitIndex[hHeart];
        const wallFit& fit = fits[w];

        const vector a
        (
            torsoShare(sigma[h], sigmaTorso, S/magS, dH, 1 - dH)
           *(S & Gi_[hHeart])
        );
        const vector e((x - heartC[hHeart])/fit.L);

        // b = J(e)^T a: grad(Vm)(x) = J(e) u for the quadratic coefficients
        // u = (a, H_xx, H_xy, H_xz, H_yy, H_yz, H_zz)
        FixedList<scalar, 9> b;
        b[0] = a.x();
        b[1] = a.y();
        b[2] = a.z();
        b[3] = a.x()*e.x();
        b[4] = a.x()*e.y() + a.y()*e.x();
        b[5] = a.x()*e.z() + a.z()*e.x();
        b[6] = a.y()*e.y();
        b[7] = a.y()*e.z() + a.z()*e.y();
        b[8] = a.z()*e.z();

        const label nRows =
            fit.valueEntries.size() + 3*fit.gradientEntries.size();
        for (label row = 0; row < nRows; ++row)
        {
            scalar cr = 0;
            for (label j = 0; j < 9; ++j)
            {
                cr += b[j]*fit.pinv(j, row);
            }
            if (row < fit.valueEntries.size())
            {
                cr *= fit.valueScale[row];
            }
            coeff.append(cr);
        }
        coeffStart.append(coeff.size());

        return w;
    };

    DynamicList<label> internalHeart, internalTorso, internalWall;
    forAll(own, faceI)
    {
        const bool ownHeart = heart[own[faceI]] > 0.5;
        if (ownHeart == (heart[nei[faceI]] > 0.5))
        {
            continue;
        }
        const label h = ownHeart ? own[faceI] : nei[faceI];
        const label t = ownHeart ? nei[faceI] : own[faceI];

        // The weight is the owner's, so the owner's distance share is
        // 1 - weight
        internalWall.append
        (
            addLayer
            (
                h,
                sigma[t],
                (ownHeart ? 1 : -1)*Sf[faceI],
                magSf[faceI],
                Cf[faceI],
                ownHeart ? 1 - weights[faceI] : weights[faceI]
            )
        );
        internalHeart.append(h);
        internalTorso.append(t);
    }

    // Coupled faces: the heart side evaluates the layer, the torso side
    // receives it
    DynamicList<label> coupledHeart, coupledSlot, coupledWall;
    DynamicList<label> receiveCell, receiveSlot;
    forAll(baseMesh_.boundary(), patchI)
    {
        const fvPatch& p = baseMesh_.boundary()[patchI];
        if (!p.coupled())
        {
            continue;
        }
        const scalarField heartN
        (
            heart.boundaryField()[patchI].patchNeighbourField()
        );
        const tensorField sigmaN
        (
            sigma.boundaryField()[patchI].patchNeighbourField()
        );
        const scalarField& pw = baseMesh_.weights().boundaryField()[patchI];
        const labelUList& cells = p.faceCells();
        const label offset = p.start() - baseMesh_.nInternalFaces();

        forAll(cells, faceI)
        {
            const label c = cells[faceI];
            if (heart[c] > 0.5 && heartN[faceI] < 0.5)
            {
                coupledWall.append
                (
                    addLayer
                    (
                        c,
                        sigmaN[faceI],
                        Sf.boundaryField()[patchI][faceI],
                        magSf.boundaryField()[patchI][faceI],
                        Cf.boundaryField()[patchI][faceI],
                        1 - pw[faceI]
                    )
                );
                coupledHeart.append(c);
                coupledSlot.append(offset + faceI);
            }
            else if (heart[c] < 0.5 && heartN[faceI] > 0.5)
            {
                receiveCell.append(c);
                receiveSlot.append(offset + faceI);
            }
        }
    }

    internalLayerHeart_.transfer(internalHeart);
    internalLayerTorso_.transfer(internalTorso);
    internalLayerWall_.transfer(internalWall);
    coupledLayerHeart_.transfer(coupledHeart);
    coupledLayerSlot_.transfer(coupledSlot);
    coupledLayerWall_.transfer(coupledWall);
    coupledReceiveCell_.transfer(receiveCell);
    coupledReceiveSlot_.transfer(receiveSlot);
    layerCoeffStart_.transfer(coeffStart);
    layerCoeff_.transfer(coeff);

    Info<< "passiveConductor: " << returnReduce(nFits, sumOp<label>())
        << " heart wall cells, quadratic fits of the wall layer built; "
        << returnReduce
           (
               internalLayerHeart_.size() + coupledLayerHeart_.size(),
               sumOp<label>()
           )
        << " heart-torso faces, "
        << returnReduce(layerCoeff_.size(), sumOp<label>())
        << " row coefficients" << endl;
}


void passiveConductor::transferWallLayer
(
    const volScalarField& Vm,
    volScalarField& rhsGlobal
) const
{
    const fvMesh& heartMesh = Vm.mesh();
    if (!wallTransferBuilt_)
    {
        buildWallTransfer(heartMesh);
        wallTransferBuilt_ = true;
    }

    // Vm and grad(Vm) over the stencil, in flat addressing
    const extendedCentredCellToCellStencil& stencil =
        centredCPCCellToCellStencilObject::New(heartMesh);
    const volVectorField gradVm(fvc::grad(Vm, "grad(Vm)"));
    const List<scalar> flatV(flatStencilField(stencil.map(), Vm));
    const List<vector> flatGrad(flatStencilField(stencil.map(), gradVm));

    // The rows of each wall cell: Vm_k - Vm_c, then grad(Vm)_k
    scalarList rows(wallRowStart_.last());
    forAll(wallCell_, w)
    {
        const scalar Vc = Vm[wallCell_[w]];
        label i = wallRowStart_[w];
        for (label k = wallValueStart_[w]; k < wallValueStart_[w + 1]; ++k)
        {
            rows[i++] = flatV[wallValueIndex_[k]] - Vc;
        }
        for
        (
            label k = wallGradientStart_[w];
            k < wallGradientStart_[w + 1];
            ++k
        )
        {
            const vector& g = flatGrad[wallGradientIndex_[k]];
            rows[i++] = g.x();
            rows[i++] = g.y();
            rows[i++] = g.z();
        }
    }

    // The layer of face f from the rows of its wall cell w
    const auto layer = [&](const label f, const label w) -> scalar
    {
        const label c0 = layerCoeffStart_[f];
        const label n = layerCoeffStart_[f + 1] - c0;
        const label r0 = wallRowStart_[w];
        scalar s = 0;
        for (label k = 0; k < n; ++k)
        {
            s += layerCoeff_[c0 + k]*rows[r0 + k];
        }
        return s;
    };

    const scalarField& V = baseMesh_.V();
    scalarField& rhs = rhsGlobal.primitiveFieldRef();

    forAll(internalLayerHeart_, i)
    {
        const scalar transfer = layer(i, internalLayerWall_[i]);
        rhs[internalLayerHeart_[i]] -= transfer/V[internalLayerHeart_[i]];
        rhs[internalLayerTorso_[i]] += transfer/V[internalLayerTorso_[i]];
    }

    // Coupled faces: the heart side evaluates the layer, the torso side
    // receives it
    scalarField boundaryLayer(baseMesh_.nBoundaryFaces(), Zero);
    const label nInternal = internalLayerHeart_.size();
    forAll(coupledLayerHeart_, i)
    {
        const scalar transfer = layer(nInternal + i, coupledLayerWall_[i]);
        rhs[coupledLayerHeart_[i]] -= transfer/V[coupledLayerHeart_[i]];
        boundaryLayer[coupledLayerSlot_[i]] = transfer;
    }
    syncTools::syncBoundaryFaceList
    (
        baseMesh_, boundaryLayer, plusEqOp<scalar>()
    );
    forAll(coupledReceiveCell_, i)
    {
        rhs[coupledReceiveCell_[i]] +=
            boundaryLayer[coupledReceiveSlot_[i]]/V[coupledReceiveCell_[i]];
    }
}


labelList passiveConductor::wallTorsoCells() const
{
    const tmp<volScalarField> tHeart
    (
        heartIndicator(baseMesh_, heartCellToBaseCell_)
    );
    const volScalarField& heart = tHeart();
    const labelUList& own = baseMesh_.owner();
    const labelUList& nei = baseMesh_.neighbour();

    boolList wallTorso(baseMesh_.nCells(), false);
    forAll(own, faceI)
    {
        const bool ownHeart = heart[own[faceI]] > 0.5;
        if (ownHeart != (heart[nei[faceI]] > 0.5))
        {
            wallTorso[ownHeart ? nei[faceI] : own[faceI]] = true;
        }
    }
    forAll(baseMesh_.boundary(), patchI)
    {
        const fvPatch& p = baseMesh_.boundary()[patchI];
        if (!p.coupled())
        {
            continue;
        }
        const scalarField heartN
        (
            heart.boundaryField()[patchI].patchNeighbourField()
        );
        forAll(p.faceCells(), faceI)
        {
            const label c = p.faceCells()[faceI];
            if (heart[c] < 0.5 && heartN[faceI] > 0.5)
            {
                wallTorso[c] = true;
            }
        }
    }

    return findIndices(wallTorso, true);
}

} // End namespace Foam

// ************************************************************************* //
