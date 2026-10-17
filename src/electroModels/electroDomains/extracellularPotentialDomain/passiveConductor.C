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
#include "PstreamReduceOps.H"
#include "fixedValueFvPatchFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"

namespace Foam
{

passiveConductor::torsoSettings passiveConductor::torsoSettings::read
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    torsoSettings s;

    s.interfaceConductivityInterpolation = dict.lookupOrDefault<word>
    (
        "interfaceConductivityInterpolation",
        "distanceWeightedHarmonic"
    );
    s.hasPhiERefPoint = dict.found("phiERefPoint");
    s.phiERefPoint =
        s.hasPhiERefPoint ? dict.get<point>("phiERefPoint") : point::zero;
    s.phiEReferenceValue =
        dict.lookupOrDefault<scalar>("phiEReferenceValue", 0.0);
    {
        const wordList zones(dict.lookup("bathCellZones"));
        s.bathCellZones = zones;
    }
    s.bathConductivityField = dict.lookupOrDefault<word>
    (
        "bathConductivityField",
        "bodyAndOrgansConductivity"
    );

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
                    << "Cannot find extracellularPotentialDomain ground patch '"
                    << s.groundPatches[patchI] << "' on mesh '" << mesh.name()
                    << "'."
                    << exit(FatalError);
            }
        }
    }

    forAll(s.surfaceCurrentPatches, patchI)
    {
        if (mesh.boundaryMesh().findPatchID(s.surfaceCurrentPatches[patchI]) < 0)
        {
            FatalErrorInFunction
                << "Cannot find extracellularPotentialDomain surface-current "
                << "patch '" << s.surfaceCurrentPatches[patchI]
                << "' on mesh '" << mesh.name() << "'."
                << exit(FatalError);
        }

        if (s.groundPatches.found(s.surfaceCurrentPatches[patchI]))
        {
            FatalErrorInFunction
                << "Extracellular-potential patch '"
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


void passiveConductor::torsoSettings::requirePhiERefPoint
(
    const word& owner
) const
{
    if (!hasPhiERefPoint)
    {
        FatalErrorInFunction
            << owner << " requires phiERefPoint only for pure "
            << "Neumann phiE boundary conditions. Add phiERefPoint, or "
            << "configure a fixedValue ground patch."
            << exit(FatalError);
    }
}


autoPtr<volScalarField> passiveConductor::readTorsoConductivity
(
    const fvMesh& mesh,
    const word& fieldName,
    const bool registerObject
)
{
    return autoPtr<volScalarField>::New
    (
        IOobject
        (
            fieldName,
            mesh.time().findInstance(mesh.dbDir(), fieldName),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            registerObject
        ),
        mesh
    );
}


passiveConductor::passiveConductor
(
    const fvMesh& baseMesh,
    const fvMesh& heartMesh,
    const labelUList& heartCellMap,
    const labelUList& heartFaceMap,
    const volTensorField& Gi,
    const volTensorField& Ge,
    const volScalarField& torsoConductivity,
    const wordList& torsoCellZones,
    const word& interfaceConductivityInterpolation,
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
    interfaceConductivityInterpolation_(interfaceConductivityInterpolation),
    sealedHeartBoundary_(sealedHeartBoundary),
    totalWallConductivity_(totalWallConductivity),
    registerFields_(registerFields),
    heartCellToBaseCell_(),
    heartFaceToBaseFace_(),
    sigmaTotalPtr_(),
    sigmaIglobalPtr_(),
    sigmaExtracellularfPtr_()
{
    if
    (
        interfaceConductivityInterpolation_ != "unweightedHarmonic"
     && interfaceConductivityInterpolation_ != "distanceWeightedHarmonic"
     && interfaceConductivityInterpolation_ != "conormalHarmonic"
    )
    {
        FatalErrorInFunction
            << "Unknown interfaceConductivityInterpolation '"
            << interfaceConductivityInterpolation_ << "'. Valid values are "
            << "unweightedHarmonic, distanceWeightedHarmonic and "
            << "conormalHarmonic."
            << exit(FatalError);
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
    assembleConductivities(Ge, torsoConductivity, torsoCellZones);
    buildSigmaExtracellularSurface();
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
    const volScalarField& torsoConductivity,
    const wordList& torsoCellZones
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

    forAll(torsoCellZones, zoneNameI)
    {
        const word& zoneName = torsoCellZones[zoneNameI];
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

} // End anonymous namespace


void passiveConductor::buildSigmaExtracellularSurface()
{
    const volTensorField& sigma = sigmaTotalPtr_();
    const volTensorField& sigmaIglobal = sigmaIglobalPtr_();
    const fvMesh& m = baseMesh_;

    sigmaExtracellularfPtr_.reset
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
                registerFields_
            ),
            m,
            dimensionedTensor("zero", sigma.dimensions(), tensor::zero)
        )
    );

    surfaceTensorField& Sef = sigmaExtracellularfPtr_();
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
        const bool ownHeart = magSqr(sigmaIntracellularI[ownCell]) > SMALL;
        const bool neiHeart = magSqr(sigmaIntracellularI[neiCell]) > SMALL;
        const bool wallFace = ownHeart != neiHeart;
        const tensor ownSigmaE =
            ownHeart && !(wallFace && totalWallConductivity_)
          ? sigmaI[ownCell] - sigmaIntracellularI[ownCell]
          : sigmaI[ownCell];
        const tensor neiSigmaE =
            neiHeart && !(wallFace && totalWallConductivity_)
          ? sigmaI[neiCell] - sigmaIntracellularI[neiCell]
          : sigmaI[neiCell];

        if (ownHeart == neiHeart)
        {
            SefI[faceI] =
                ownHeart
              ? extracellularFaceConductivity::linear
                (
                    ownSigmaE, neiSigmaE, weights[faceI]
                )
              : torsoFaceConductivity
                (
                    ownSigmaE, neiSigmaE, weights[faceI]
                );
        }
        else if (interfaceConductivityInterpolation_ == "unweightedHarmonic")
        {
            SefI[faceI] = extracellularFaceConductivity::unweightedHarmonic
            (
                ownSigmaE, neiSigmaE
            );
        }
        else if (interfaceConductivityInterpolation_ == "conormalHarmonic")
        {
            SefI[faceI] = extracellularFaceConductivity::conormalHarmonic
            (
                ownSigmaE,
                neiSigmaE,
                m.Sf()[faceI]/m.magSf()[faceI],
                1.0 - weights[faceI],
                weights[faceI]
            );
        }
        else
        {
            SefI[faceI] =
                extracellularFaceConductivity::distanceWeightedHarmonic
                (
                    ownSigmaE,
                    neiSigmaE,
                    1.0 - weights[faceI],
                    weights[faceI]
                );
        }
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

            Field<tensor>& Sefp = Sef.boundaryFieldRef()[patchI];
            const scalarField& patchWeights =
                m.weights().boundaryField()[patchI];
            const vectorField patchNf(m.boundary()[patchI].nf());

            forAll(Sefp, faceI)
            {
                const bool pHeart = magSqr(sigmaIP[faceI]) > SMALL;
                const bool nHeart = magSqr(sigmaIN[faceI]) > SMALL;
                const bool wallFace = pHeart != nHeart;
                const tensor pSigmaE =
                    pHeart && !(wallFace && totalWallConductivity_)
                  ? sigmaP[faceI] - sigmaIP[faceI]
                  : sigmaP[faceI];
                const tensor nSigmaE =
                    nHeart && !(wallFace && totalWallConductivity_)
                  ? sigmaN[faceI] - sigmaIN[faceI]
                  : sigmaN[faceI];
                if (pHeart == nHeart)
                {
                    Sefp[faceI] =
                        pHeart
                      ? extracellularFaceConductivity::linear
                        (
                            pSigmaE, nSigmaE, patchWeights[faceI]
                        )
                      : torsoFaceConductivity
                        (
                            pSigmaE, nSigmaE, patchWeights[faceI]
                        );
                }
                else if
                (
                    interfaceConductivityInterpolation_
                 == "unweightedHarmonic"
                )
                {
                    Sefp[faceI] =
                        extracellularFaceConductivity::unweightedHarmonic
                        (
                            pSigmaE, nSigmaE
                        );
                }
                else if
                (
                    interfaceConductivityInterpolation_
                 == "conormalHarmonic"
                )
                {
                    Sefp[faceI] =
                        extracellularFaceConductivity::conormalHarmonic
                        (
                            pSigmaE,
                            nSigmaE,
                            patchNf[faceI],
                            1.0 - patchWeights[faceI],
                            patchWeights[faceI]
                        );
                }
                else
                {
                    Sefp[faceI] = extracellularFaceConductivity::
                        distanceWeightedHarmonic
                        (
                            pSigmaE,
                            nSigmaE,
                            1.0 - patchWeights[faceI],
                            patchWeights[faceI]
                        );
                }
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
}


tmp<fvScalarMatrix> passiveConductor::assemble
(
    volScalarField& psi,
    volScalarField& heartPsi,
    const volScalarField* source,
    faceCorrection* correction
) const
{
    tmp<fvScalarMatrix> tphiEqn;
    if (interfaceConductivityInterpolation_ == "conormalHarmonic")
    {
        const surfaceTensorField& Sef = sigmaExtracellularfPtr_();
        const surfaceVectorField& Sf = baseMesh_.Sf();
        const surfaceScalarField gammaNN
        (
            (Sf & Sef & Sf)/sqr(baseMesh_.magSf())
        );
        const surfaceVectorField SfGammaCorr((Sf & Sef) - gammaNN*Sf);
        // Heart-side wall conductivity for the two-sided flux continuity
        const volTensorField sigmaWall
        (
            totalWallConductivity_
          ? volTensorField(sigmaTotalPtr_())
          : volTensorField(sigmaTotalPtr_() - sigmaIglobalPtr_())
        );
        tmp<volVectorField> tgrad(fvc::grad(psi));
        for (label iter = 0; iter < 2; ++iter)
        {
            tgrad = interfaceGradient
            (
                psi, sigmaWall, sigmaIglobalPtr_(), tgrad()
            );
        }
        tphiEqn =
        (
            fvm::laplacian(gammaNN, psi)
          + fvc::div(SfGammaCorr & linearInterpolate(tgrad()))
        );
        if (source)
        {
            tphiEqn.ref() -= *source;
        }

        // Non-orthogonal correction of fvm::laplacian rebuilt from the jump-aware gradient
        const surfaceScalarField nonOrthGradCorr
        (
            gammaNN*baseMesh_.magSf()
           *(
                baseMesh_.nonOrthCorrectionVectors()
              & linearInterpolate(tgrad() - fvc::grad(psi))
            )
        );
        tphiEqn.ref() += fvc::div(nonOrthGradCorr);

        if (correction)
        {
            correction->crossFlux = SfGammaCorr & linearInterpolate(tgrad());
            correction->crossFlux.ref() += nonOrthGradCorr;
            correction->gradient = tgrad;
        }
    }
    else
    {
        tphiEqn = fvm::laplacian(sigmaExtracellularfPtr_(), psi);
        if (source)
        {
            tphiEqn.ref() -= *source;
        }

        if (correction)
        {
            correction->gradient = fvc::grad(psi);
        }
    }
    fvScalarMatrix& phiEqn = tphiEqn.ref();

    {
        fvScalarMatrix heartPhiEqn
        (
            fvm::laplacian
            (
                insulatedFaceConductivity
                (
                    Gi_, sealingField_, sealedHeartBoundary_
                ),
                heartPsi
            )
        );
        const labelUList& heartFaceMap = heartFaceToBaseFace_;
        const labelUList& heartOwner = heartPsi.mesh().owner();

        forAll(heartCellToBaseCell_, heartCellI)
        {
            const label baseCellI = heartCellToBaseCell_[heartCellI];
            phiEqn.diag()[baseCellI] += heartPhiEqn.diag()[heartCellI];
            phiEqn.source()[baseCellI] += heartPhiEqn.source()[heartCellI];
        }
        forAll(heartOwner, heartFaceI)
        {
            const label baseFaceI = heartFaceMap[heartFaceI];
            phiEqn.upper()[baseFaceI] += heartPhiEqn.upper()[heartFaceI];
        }

        forAll(heartPhiEqn.internalCoeffs(), patchI)
        {
            const fvPatch& heartPatch = heartPsi.mesh().boundary()[patchI];
            forAll(heartPhiEqn.internalCoeffs()[patchI], faceI)
            {
                const scalar internalCoeff =
                    heartPhiEqn.internalCoeffs()[patchI][faceI];
                const scalar boundaryCoeff =
                    heartPhiEqn.boundaryCoeffs()[patchI][faceI];
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
                phiEqn.internalCoeffs()[basePatchI][basePatchFaceI] +=
                    internalCoeff;
                phiEqn.boundaryCoeffs()[basePatchI][basePatchFaceI] +=
                    boundaryCoeff;
            }
        }
    }

    return tphiEqn;
}


void passiveConductor::applyReference
(
    fvScalarMatrix& eqn,
    const volScalarField& psi,
    const label refCell,
    const scalar value
) const
{
    // Only pin a reference cell when the BC inventory leaves the system
    // singular (pure Neumann). When any psi patch is fixedValue the
    // Dirichlet patch already determines the constant; an extra pin would
    // over-determine the discrete problem with a value that need not match
    // the analytical psi at that cell, polluting the entire field by
    // O(ref offset).
    const auto& bf = psi.boundaryField();
    forAll(bf, patchI)
    {
        if (isA<fixedValueFvPatchScalarField>(bf[patchI]))
        {
            return;
        }
    }

    if (refCell >= 0)
    {
        eqn.setReference(refCell, value, true);
    }
}


label passiveConductor::referenceCell(const point& p) const
{
    const label refCell = baseMesh_.findCell(p);
    const label localOwns = refCell >= 0 ? 1 : 0;
    label ownerCount = localOwns;
    reduce(ownerCount, sumOp<label>());

    if (ownerCount != 1)
    {
        FatalErrorInFunction
            << "phiERefPoint " << p
            << " must be owned by exactly one mesh partition, but "
            << ownerCount << " partitions reported a containing cell."
            << exit(FatalError);
    }

    return refCell;
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
