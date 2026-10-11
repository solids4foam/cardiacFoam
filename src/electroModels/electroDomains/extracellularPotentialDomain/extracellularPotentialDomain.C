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

#include "extracellularPotentialDomain.H"
#include "extracellularFaceConductivity.H"
#include "myocardiumDomainInterface.H"
#include "electroVerificationModel.H"
#include "insulatedFaceConductivity.H"
#include "sealedWallPolicy.H"
#include "exposedPhiETrace.H"
#include "interfaceGradient.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"
#include "OStringStream.H"
#include "dimVoltage.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "nonOrthogonalCorrectorLoop.H"

namespace Foam
{

defineTypeNameAndDebug(extracellularPotentialDomain, 0);
addToRunTimeSelectionTable
(
    electroStateDomain,
    extracellularPotentialDomain,
    dictionary
);

namespace
{

// Resolve the non-orthogonal corrector count from the PIMPLE dictionary.
label resolveNonOrthogonalCorrectors(const fvMesh& baseMesh)
{
    const dictionary& pimpleDict =
        baseMesh.solutionDict().subOrEmptyDict("PIMPLE");

    const label nCorr =
        pimpleDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);

    if (nCorr < 0)
    {
        FatalIOErrorInFunction(pimpleDict)
            << "nNonOrthogonalCorrectors must be non-negative; found "
            << nCorr << '.' << exit(FatalIOError);
    }

    return nCorr;
}


wordList potentialPatchTypes(const fvMesh& mesh, const dictionary& dict)
{
    wordList patchTypes
    (
        mesh.boundary().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    if (const dictionary* groundDict = dict.findDict("groundPatches"))
    {
        const wordList patchNames(groundDict->toc());

        forAll(patchNames, patchI)
        {
            const label patchId =
                mesh.boundaryMesh().findPatchID(patchNames[patchI]);

            if (patchId < 0)
            {
                FatalErrorInFunction
                    << "Cannot find extracellularPotentialDomain ground patch '"
                    << patchNames[patchI] << "' on mesh '" << mesh.name()
                    << "'."
                    << exit(FatalError);
            }

            patchTypes[patchId] = fixedValueFvPatchScalarField::typeName;
        }
    }

    if (const dictionary* currentDict = dict.findDict("surfaceCurrentPatches"))
    {
        const wordList patchNames(currentDict->toc());

        forAll(patchNames, patchI)
        {
            const label patchId =
                mesh.boundaryMesh().findPatchID(patchNames[patchI]);

            if (patchId < 0)
            {
                FatalErrorInFunction
                    << "Cannot find extracellularPotentialDomain surface-current "
                    << "patch '" << patchNames[patchI] << "' on mesh '"
                    << mesh.name() << "'."
                    << exit(FatalError);
            }

            if
            (
                patchTypes[patchId]
             == fixedValueFvPatchScalarField::typeName
            )
            {
                FatalErrorInFunction
                    << "Extracellular-potential patch '" << patchNames[patchI]
                    << "' cannot be listed in both groundPatches and "
                    << "surfaceCurrentPatches."
                    << exit(FatalError);
            }

            patchTypes[patchId] = fixedGradientFvPatchScalarField::typeName;
        }
    }

    return patchTypes;
}


void applyGroundPatchValues(volScalarField& phiE, const dictionary& dict)
{
    const dictionary* groundDict = dict.findDict("groundPatches");

    if (!groundDict)
    {
        return;
    }

    const wordList patchNames(groundDict->toc());

    forAll(patchNames, patchI)
    {
        const word& patchName = patchNames[patchI];
        const label patchId = phiE.mesh().boundaryMesh().findPatchID(patchName);

        const scalar value = groundDict->get<scalar>(patchName);
        phiE.boundaryFieldRef()[patchId] = value;
    }
}


} // End anonymous namespace


extracellularPotentialDomain::extracellularPotentialDomain
(
    const fvMesh& baseMesh,
    myocardiumDomainInterface& heartDomain,
    const dictionary& dict
)
:
    baseMesh_(baseMesh),
    heartDomain_(heartDomain),
    sigmaTotalPtr_(),
    sigmaIglobalPtr_(),
    sigmaExtracellularfPtr_(),
    interfaceConductivityInterpolation_
    (
        dict.lookupOrDefault<word>
        (
            "interfaceConductivityInterpolation",
            "distanceWeightedHarmonic"
        )
    ),
    phiEPtr_(),
    VmGlobalPtr_(),
    heartCellToBaseCell_(),
    phiEReferencePoint_
    (
        dict.found("phiERefPoint")
      ? dict.get<point>("phiERefPoint")
      : point::zero
    ),
    phiEReferenceValue_(dict.lookupOrDefault<scalar>("phiEReferenceValue", 0.0)),
    hasPhiEReferencePoint_(dict.found("phiERefPoint")),
    bathCellZoneNames_(dict.lookup("bathCellZones")),
    bathConductivityFieldName_
    (
        dict.lookupOrDefault<word>
        (
            "bathConductivityField",
            "bodyAndOrgansConductivity"
        )
    ),
    surfaceCurrentPatchNames_(),
    surfaceCurrentPatchValues_(),
    hasDirichletPatch_(false),
    nNonOrthogonalCorrectors_
    (
        resolveNonOrthogonalCorrectors(baseMesh)
    ),
    sealedHeartBoundary_
    (
        dict.parent().get<Switch>("sealedHeartBoundary")
    ),
    manufacturedEvaluationTime_(baseMesh.time().value())
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
    if (const dictionary* currentDict = dict.findDict("surfaceCurrentPatches"))
    {
        surfaceCurrentPatchNames_ = currentDict->toc();
        surfaceCurrentPatchValues_.setSize
        (
            surfaceCurrentPatchNames_.size(),
            0.0
        );

        forAll(surfaceCurrentPatchNames_, patchI)
        {
            surfaceCurrentPatchValues_[patchI] =
                currentDict->get<scalar>(surfaceCurrentPatchNames_[patchI]);
        }
    }

    phiEPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "phiE",
                baseMesh_.time().timeName(),
                baseMesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            baseMesh_,
            dimensionedScalar("phiE", dimVoltage, 0.0),
            potentialPatchTypes(baseMesh_, dict)
        )
    );
    {
        const wordList expectedTypes(potentialPatchTypes(baseMesh_, dict));
        const auto& bf = phiEPtr_().boundaryField();
        OStringStream mismatches;
        label nMismatched = 0;

        forAll(bf, patchI)
        {
            const polyPatch& pp = baseMesh_.boundaryMesh()[patchI];
            if (polyPatch::constraintType(pp.type()))
            {
                continue;
            }
            if (bf[patchI].type() != expectedTypes[patchI])
            {
                mismatches
                    << ' ' << pp.name() << " (" << bf[patchI].type()
                    << ", expected " << expectedTypes[patchI] << ')';
                ++nMismatched;
            }
        }

        if (nMismatched)
        {
            FatalErrorInFunction
                << "The global phiE patch types differ from those set by "
                << "groundPatches, surfaceCurrentPatches and the zeroGradient "
                << "default:" << mismatches.str().c_str()
                << ". A phiE file in the start time directory sets them; "
                << "remove it or give these patches the expected types."
                << exit(FatalError);
        }
    }

    applyGroundPatchValues(phiEPtr_(), dict);
    sealedWallPolicy::publish(phiEPtr_(), false);

    {
        const auto& bf = phiEPtr_().boundaryField();
        forAll(bf, patchI)
        {
            if (isA<fixedValueFvPatchScalarField>(bf[patchI]))
            {
                hasDirichletPatch_ = true;
                break;
            }
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
                IOobject::NO_WRITE
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
                IOobject::NO_WRITE
            ),
            baseMesh_,
            dimensionedTensor("zero", conductivityDim, tensor::zero),
            "zeroGradient"
        )
    );

    VmGlobalPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "VmGlobal",
                baseMesh_.time().timeName(),
                baseMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            baseMesh_,
            dimensionedScalar("VmGlobal", dimVoltage, 0.0),
            "zeroGradient"
        )
    );

    if (const electroVerificationModel* model = verifier())
    {
        const wordList claimed(model->boundaryDataPatches());
        if (!claimed.empty())
        {
            wordHashSet settable;
            if (const dictionary* d = dict.findDict("groundPatches"))
            {
                settable.insert(d->toc());
            }
            if (const dictionary* d = dict.findDict("surfaceCurrentPatches"))
            {
                settable.insert(d->toc());
            }

            const wordHashSet claimedSet(claimed);
            if (claimedSet.size() != claimed.size())
            {
                FatalErrorInFunction
                    << "Verification boundary-data patches " << claimed
                    << " contain a duplicate name."
                    << exit(FatalError);
            }
            if (claimedSet != settable)
            {
                FatalErrorInFunction
                    << "Verification boundary-data patches " << claimed
                    << " must equal the groundPatches and "
                    << "surfaceCurrentPatches entries "
                    << settable.sortedToc() << '.'
                    << exit(FatalError);
            }
        }

        if
        (
            hasDirichletPatch_
         && model->requiresPhiECompatibilityProjection()
        )
        {
            WarningInFunction
                << "The phiE compatibility projection is ignored because a "
                << "fixedValue patch determines the constant."
                << endl;
        }

        if (model->requiresInterfaceFluxDiagnostic())
        {
            if (!baseMesh_.fluxRequired(phiEPtr_().name()))
            {
                FatalErrorInFunction
                    << "The interface-flux diagnostic needs phiE in the "
                    << "fluxRequired sub-dictionary of fvSchemes."
                    << exit(FatalError);
            }

            interfaceFluxPtr_.reset
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "phiEInterfaceFlux",
                        baseMesh_.time().timeName(),
                        baseMesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    baseMesh_,
                    dimensionedScalar(dimCurrent, Zero)
                )
            );
            interfaceFluxPtr_().setOriented();

            interfaceGradientPtr_.reset
            (
                new volVectorField
                (
                    IOobject
                    (
                        "phiEInterfaceGradient",
                        baseMesh_.time().timeName(),
                        baseMesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    baseMesh_,
                    dimensionedVector
                    (
                        phiEPtr_().dimensions()/dimLength, Zero
                    )
                )
            );
        }
    }

    buildHeartScatterMap();
    assembleConductivities();
    buildSigmaExtracellularSurface();
    heartDomain_.bindExternalPhiE(phiEPtr_(), heartCellToBaseCell_);
}

void extracellularPotentialDomain::buildHeartScatterMap()
{
    const labelUList* heartCellMapPtr = heartDomain_.subsetCellMapPtr();
    const fvMesh& heartMesh =
        static_cast<const electroStateProvider&>(heartDomain_).mesh();

    if (!heartCellMapPtr)
    {
        FatalErrorInFunction
            << "Unified phiE requires the myocardium domain to expose a "
            << "subset cell map to the base mesh."
            << exit(FatalError);
    }

    heartCellToBaseCell_ = labelList(*heartCellMapPtr);

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
        const labelUList* heartFaceMapPtr = heartDomain_.subsetFaceMapPtr();

        if (!heartFaceMapPtr)
        {
            FatalErrorInFunction
                << "Matched intracellular assembly requires the heart/base "
                << "face map."
                << exit(FatalError);
        }

        const labelUList& heartFaceMap = *heartFaceMapPtr;
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
    }
}

void extracellularPotentialDomain::assembleConductivities()
{
    const volTensorField* GiPtr = heartDomain_.intracellularConductivityPtr();
    const volTensorField* GePtr = heartDomain_.extracellularConductivityPtr();

    if (!GiPtr || !GePtr)
    {
        FatalErrorInFunction
            << "extracellularPotentialDomain requires split intracellular and "
            << "extracellular conductivities from a bidomain myocardium."
            << exit(FatalError);
    }

    volTensorField& sigmaTotal = sigmaTotalPtr_();
    volTensorField& sigmaIglobal = sigmaIglobalPtr_();
    tensorField& sigmaTotalI = sigmaTotal.primitiveFieldRef();
    tensorField& sigmaIglobalI = sigmaIglobal.primitiveFieldRef();
    const tensorField& Gi = GiPtr->primitiveField();
    const tensorField& Ge = GePtr->primitiveField();

    if
    (
        heartCellToBaseCell_.size() != Gi.size()
     || heartCellToBaseCell_.size() != Ge.size()
    )
    {
        FatalErrorInFunction
            << "Heart/base cell map size " << heartCellToBaseCell_.size()
            << " does not match intracellular/extracellular conductivity "
            << "cell counts " << Gi.size() << "/" << Ge.size() << "."
            << exit(FatalError);
    }

    labelList materialRegion(baseMesh_.nCells(), 0);

    forAll(heartCellToBaseCell_, heartCellI)
    {
        const label baseCellI = heartCellToBaseCell_[heartCellI];
        const tensor heartSigma = Gi[heartCellI] + Ge[heartCellI];

        if (magSqr(heartSigma) <= SMALL)
        {
            FatalErrorInFunction
                << "Myocardium cell " << heartCellI
                << " has zero total conductivity."
                << exit(FatalError);
        }

        sigmaTotalI[baseCellI] = heartSigma;
        sigmaIglobalI[baseCellI] = Gi[heartCellI];
        materialRegion[baseCellI] = 1;
    }

    volScalarField bathSigma
    (
        IOobject
        (
            bathConductivityFieldName_,
            baseMesh_.time().findInstance
            (
                baseMesh_.dbDir(),
                bathConductivityFieldName_
            ),
            baseMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        baseMesh_
    );

    const scalarField& bathSigmaI = bathSigma.primitiveField();

    forAll(bathCellZoneNames_, zoneNameI)
    {
        const word& zoneName = bathCellZoneNames_[zoneNameI];
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

void extracellularPotentialDomain::buildSigmaExtracellularSurface()
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
                IOobject::NO_WRITE
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
        const tensor ownSigmaE =
            ownHeart
          ? sigmaI[ownCell] - sigmaIntracellularI[ownCell]
          : sigmaI[ownCell];
        const tensor neiSigmaE =
            neiHeart
          ? sigmaI[neiCell] - sigmaIntracellularI[neiCell]
          : sigmaI[neiCell];

        if (ownHeart == neiHeart)
        {
            SefI[faceI] = extracellularFaceConductivity::linear
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
                const tensor pSigmaE = pHeart
                  ? sigmaP[faceI] - sigmaIP[faceI]
                  : sigmaP[faceI];
                const tensor nSigmaE = nHeart
                  ? sigmaN[faceI] - sigmaIN[faceI]
                  : sigmaN[faceI];
                if (pHeart == nHeart)
                {
                    Sefp[faceI] = extracellularFaceConductivity::linear
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
            // Physical boundaries use the cell-internal extracellular value.
            Sef.boundaryFieldRef()[patchI] =
                sigma.boundaryField()[patchI].patchInternalField()
              - sigmaIglobal.boundaryField()[patchI].patchInternalField();
        }
    }
}


label extracellularPotentialDomain::referenceCell() const
{
    if (!hasPhiEReferencePoint_)
    {
        FatalErrorInFunction
            << "extracellularPotentialDomain requires phiERefPoint only for pure "
            << "Neumann phiE boundary conditions. Add phiERefPoint, or "
            << "configure a fixedValue ground patch."
            << exit(FatalError);
    }

    const label refCell = baseMesh_.findCell(phiEReferencePoint_);
    const label localOwns = refCell >= 0 ? 1 : 0;
    label ownerCount = localOwns;
    reduce(ownerCount, sumOp<label>());

    if (ownerCount != 1)
    {
        FatalErrorInFunction
            << "phiERefPoint " << phiEReferencePoint_
            << " must be owned by exactly one mesh partition, but "
            << ownerCount << " partitions reported a containing cell."
            << exit(FatalError);
    }

    return refCell;
}

void extracellularPotentialDomain::scatterHeartVm()
{
    const volScalarField* heartVmPtr = heartDomain_.VmPtr();

    if (!heartVmPtr)
    {
        FatalErrorInFunction
            << "Heart domain does not expose Vm."
            << exit(FatalError);
    }

    scalarField& VmGlobal = VmGlobalPtr_().primitiveFieldRef();
    const scalarField& heartVm = heartVmPtr->primitiveField();
    VmGlobal = 0.0;

    forAll(heartCellToBaseCell_, heartCellI)
    {
        VmGlobal[heartCellToBaseCell_[heartCellI]] = heartVm[heartCellI];
    }

    VmGlobalPtr_().correctBoundaryConditions();
}

const electroVerificationModel* extracellularPotentialDomain::verifier() const
{
    return heartDomain_.verificationModel();
}


void extracellularPotentialDomain::applyManufacturedBoundaryData()
{
    const electroVerificationModel* model = verifier();
    if (!model)
    {
        return;
    }

    const wordList patchNames(model->boundaryDataPatches());
    if (patchNames.empty())
    {
        return;
    }

    volScalarField& phiE = phiEPtr_();
    model->applyPhiEBoundaryData(phiE, manufacturedEvaluationTime_);

    for (const word& patchName : patchNames)
    {
        phiE.boundaryFieldRef()
        [
            baseMesh_.boundaryMesh().findPatchID(patchName)
        ].evaluate();
    }
}


void extracellularPotentialDomain::projectPhiECompatibility
(
    fvScalarMatrix& phiEqn
) const
{
    scalar imbalance = sum(phiEqn.source());
    forAll(phiEqn.boundaryCoeffs(), patchI)
    {
        if (!baseMesh_.boundary()[patchI].coupled())
        {
            imbalance += sum(phiEqn.boundaryCoeffs()[patchI]);
        }
    }
    reduce(imbalance, sumOp<scalar>());

    const scalarField& V = baseMesh_.V();
    phiEqn.source() -= V*(imbalance/gSum(V));

    Info<< "phiE compatibility projection: imbalance " << imbalance << nl;
}


void extracellularPotentialDomain::storeInterfaceFlux
(
    const surfaceScalarField& F
)
{
    surfaceScalarField& flux = interfaceFluxPtr_();
    flux = dimensionedScalar(flux.dimensions(), Zero);

    const volTensorField& sigmaI = sigmaIglobalPtr_();
    const labelUList& owner = baseMesh_.owner();
    const labelUList& neighbour = baseMesh_.neighbour();
    auto isHeart = [](const tensor& g) { return magSqr(g) > SMALL; };

    forAll(owner, faceI)
    {
        if
        (
            isHeart(sigmaI[owner[faceI]])
         != isHeart(sigmaI[neighbour[faceI]])
        )
        {
            flux[faceI] = F[faceI];
        }
    }

    forAll(flux.boundaryField(), patchI)
    {
        if (!baseMesh_.boundary()[patchI].coupled())
        {
            continue;
        }

        const tensorField sigmaIP
        (
            sigmaI.boundaryField()[patchI].patchInternalField()
        );
        const tensorField sigmaIN
        (
            sigmaI.boundaryField()[patchI].patchNeighbourField()
        );
        const scalarField& Fp = F.boundaryField()[patchI];
        scalarField& pFlux = flux.boundaryFieldRef()[patchI];

        forAll(pFlux, faceI)
        {
            if (isHeart(sigmaIP[faceI]) != isHeart(sigmaIN[faceI]))
            {
                pFlux[faceI] = Fp[faceI];
            }
        }
    }
}


void extracellularPotentialDomain::prepareTimeStep(scalar t0, scalar dt)
{
    manufacturedEvaluationTime_ = t0 + dt;
    applyManufacturedBoundaryData();
}


void extracellularPotentialDomain::solvePhiEOnce()
{
    forAll(surfaceCurrentPatchNames_, patchI)
    {
        const label patchId =
            baseMesh_.boundaryMesh().findPatchID(surfaceCurrentPatchNames_[patchI]);

        if (patchId < 0)
        {
            FatalErrorInFunction
                << "Cannot find extracellularPotentialDomain surface-current patch '"
                << surfaceCurrentPatchNames_[patchI] << "' on mesh '"
                << baseMesh_.name() << "'."
                << exit(FatalError);
        }

        fixedGradientFvPatchScalarField& phiEPatch =
            refCast<fixedGradientFvPatchScalarField>
            (
                phiEPtr_().boundaryFieldRef()[patchId]
            );

        tmp<tensorField> tSigmaPatch =
            sigmaTotalPtr_().boundaryField()[patchId].patchInternalField();
        const tensorField& sigmaPatch = tSigmaPatch();
        tmp<vectorField> tNf(phiEPatch.patch().nf());
        const vectorField& nf = tNf();
        scalarField& gradient = phiEPatch.gradient();

        forAll(gradient, faceI)
        {
            const scalar sigmaN =
                nf[faceI] & (sigmaPatch[faceI] & nf[faceI]);

            if (sigmaN <= SMALL)
            {
                FatalErrorInFunction
                    << "Surface-current patch '"
                    << surfaceCurrentPatchNames_[patchI]
                    << "' found non-positive normal conductivity " << sigmaN
                    << " on face " << faceI << "."
                    << exit(FatalError);
            }

            gradient[faceI] = surfaceCurrentPatchValues_[patchI]/sigmaN;
        }

        phiEPatch.evaluate();
    }

    applyManufacturedBoundaryData();

    const volTensorField* GiPtr = heartDomain_.intracellularConductivityPtr();
    const volScalarField* VmPtr = heartDomain_.VmPtr();

    if (!GiPtr || !VmPtr)
    {
        FatalErrorInFunction
            << "extracellularPotentialDomain requires heart Vm and intracellular "
            << "conductivity fields."
            << exit(FatalError);
    }

    volScalarField* heartPhiEPtr = const_cast<volScalarField*>
    (
        heartDomain_.phiEPtr()
    );
    if (!heartPhiEPtr)
    {
        FatalErrorInFunction
            << "Matched intracellular assembly requires the heart phiE field."
            << exit(FatalError);
    }

    setExposedPhiETrace
    (
        *heartPhiEPtr, phiEPtr_(), heartCellToBaseCell_,
        *heartDomain_.subsetFaceMapPtr(),
        *heartDomain_.extracellularConductivityPtr(),
        sigmaTotalPtr_()
    );

    tmp<volScalarField> tHeartRhs =
        -fvc::laplacian
        (
            insulatedFaceConductivity(*GiPtr, *VmPtr, sealedHeartBoundary_),
            *VmPtr
        );
    const volScalarField& heartRhs = tHeartRhs();

    volScalarField rhsGlobal
    (
        IOobject
        (
            "phiESource",
            baseMesh_.time().timeName(),
            baseMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        baseMesh_,
        dimensionedScalar("zero", heartRhs.dimensions(), 0.0),
        "zeroGradient"
    );

    scalarField& rhsGlobalI = rhsGlobal.primitiveFieldRef();
    const scalarField& heartRhsI = heartRhs.primitiveField();
    rhsGlobalI = 0.0;

    forAll(heartCellToBaseCell_, heartCellI)
    {
        rhsGlobalI[heartCellToBaseCell_[heartCellI]] = heartRhsI[heartCellI];
    }

    rhsGlobal.correctBoundaryConditions();

    if (const electroVerificationModel* model = verifier())
    {
        model->addManufacturedPhiESource
        (
            rhsGlobal,
            manufacturedEvaluationTime_
        );
    }

    const bool storeFlux =
        verifier() && verifier()->requiresInterfaceFluxDiagnostic();
    tmp<surfaceScalarField> tCrossFlux;

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
        const volTensorField sigmaE(sigmaTotalPtr_() - sigmaIglobalPtr_());
        tmp<volVectorField> tgrad(fvc::grad(phiEPtr_()));
        for (label iter = 0; iter < 2; ++iter)
        {
            tgrad = interfaceGradient
            (
                phiEPtr_(), sigmaE, sigmaIglobalPtr_(), tgrad()
            );
        }
        tphiEqn =
        (
            fvm::laplacian(gammaNN, phiEPtr_())
          + fvc::div(SfGammaCorr & linearInterpolate(tgrad()))
         == rhsGlobal
        );

        // Non-orthogonal correction of fvm::laplacian rebuilt from the jump-aware gradient
        const surfaceScalarField nonOrthGradCorr
        (
            gammaNN*baseMesh_.magSf()
           *(
                baseMesh_.nonOrthCorrectionVectors()
              & linearInterpolate(tgrad() - fvc::grad(phiEPtr_()))
            )
        );
        tphiEqn.ref() += fvc::div(nonOrthGradCorr);

        if (storeFlux)
        {
            tCrossFlux = SfGammaCorr & linearInterpolate(tgrad());
            tCrossFlux.ref() += nonOrthGradCorr;
            interfaceGradientPtr_() = tgrad();
        }
    }
    else
    {
        tphiEqn =
        (
            fvm::laplacian(sigmaExtracellularfPtr_(), phiEPtr_())
         == rhsGlobal
        );

        if (storeFlux)
        {
            interfaceGradientPtr_() = fvc::grad(phiEPtr_());
        }
    }
    fvScalarMatrix& phiEqn = tphiEqn.ref();

    {
        fvScalarMatrix heartPhiEqn
        (
            fvm::laplacian
            (
                insulatedFaceConductivity(*GiPtr, *VmPtr, sealedHeartBoundary_),
                *heartPhiEPtr
            )
        );
        const labelUList* heartFaceMapPtr = heartDomain_.subsetFaceMapPtr();
        if (!heartFaceMapPtr)
        {
            FatalErrorInFunction
                << "Matched intracellular assembly requires the heart/base "
                << "face map."
                << exit(FatalError);
        }
        const labelUList& heartFaceMap = *heartFaceMapPtr;
        const labelUList& heartOwner = heartPhiEPtr->mesh().owner();

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
            const fvPatch& heartPatch = heartPhiEPtr->mesh().boundary()[patchI];
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

    if
    (
        !hasDirichletPatch_
     && verifier()
     && verifier()->requiresPhiECompatibilityProjection()
    )
    {
        projectPhiECompatibility(phiEqn);
    }

    // Only pin a reference cell when the BC inventory leaves the system
    // singular (pure Neumann). When any phiE patch is fixedValue the
    // Dirichlet patch already determines the constant; an extra pin would
    // over-determine the discrete problem with a value (phiEReferenceValue_)
    // that need not match the analytical phiE at that cell, polluting the
    // entire field by O(ref offset).
    if (!hasDirichletPatch_)
    {
        const label refCell = referenceCell();
        if (refCell >= 0)
        {
            phiEqn.setReference(refCell, phiEReferenceValue_, true);
        }
    }

    solve(phiEqn);

    if (storeFlux)
    {
        tmp<surfaceScalarField> tF(phiEqn.flux());
        if (tCrossFlux)
        {
            tF.ref() += tCrossFlux();
        }
        storeInterfaceFlux(tF());
    }
}


void extracellularPotentialDomain::solvePhiE()
{
    scatterHeartVm();

    // Reassemble and solve global phiE with the configured correction count.
    correctNonOrthogonalLoop
    (
        nNonOrthogonalCorrectors_,
        [&]() { solvePhiEOnce(); }
    );
}

void extracellularPotentialDomain::advance(scalar t0, scalar dt)
{
    (void)t0; (void)dt;
    solvePhiE();
}

void extracellularPotentialDomain::write()
{
    if (phiEPtr_.valid()) phiEPtr_().write();
    if (sigmaTotalPtr_.valid()) sigmaTotalPtr_().write();
    if (VmGlobalPtr_.valid()) VmGlobalPtr_().write();
}

} // End namespace Foam

// ************************************************************************* //
