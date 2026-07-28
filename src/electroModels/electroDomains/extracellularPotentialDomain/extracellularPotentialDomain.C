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
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"
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

// Resolve the phiE non-orthogonal corrector count. Default: the same
// nNonOrthogonalCorrectors Vm uses (system/fvSolution PIMPLE), so forgetting
// it can never silently drop phiE to zero correction. Optional co-located
// override: PIMPLE/phiENonOrthogonalCorrectors. Fail loud if the key is left
// in the retired constant/electroProperties -> bathPotentialDomain location.
label resolvePhiENonOrthogonalCorrectors
(
    const fvMesh& baseMesh,
    const dictionary& dict
)
{
    if (dict.found("phiENonOrthogonalCorrectors"))
    {
        FatalIOErrorInFunction(dict)
            << "phiENonOrthogonalCorrectors has moved to system/fvSolution"
            << " (PIMPLE), co-located with nNonOrthogonalCorrectors." << nl
            << "Remove it from constant/electroProperties -> bathPotentialDomain"
            << " and, only if you need phiE to differ from Vm, set it under"
            << " PIMPLE." << exit(FatalIOError);
    }

    const dictionary& pimpleDict =
        baseMesh.solutionDict().subOrEmptyDict("PIMPLE");
    const label defaultNCorr =
        pimpleDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);
    const label nCorr =
        pimpleDict.lookupOrDefault<label>
        (
            "phiENonOrthogonalCorrectors",
            defaultNCorr
        );

    if (nCorr < 0)
    {
        FatalIOErrorInFunction(pimpleDict)
            << "phiENonOrthogonalCorrectors must be non-negative; found "
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
    sigmaTotalfPtr_(),
    sigmaExtracellularfPtr_(),
    interfaceConductivityInterpolation_
    (
        dict.lookupOrDefault<word>
        (
            "interfaceConductivityInterpolation",
            "distanceWeightedHarmonic"
        )
    ),
    intracellularAssembly_
    (
        dict.lookupOrDefault<word>("intracellularAssembly", "matchedSubmesh")
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
    reportSetup_(dict.lookupOrDefault<Switch>("reportSetup", false)),
    phiENonOrthogonalCorrectors_
    (
        resolvePhiENonOrthogonalCorrectors(baseMesh, dict)
    )
{
    if
    (
        interfaceConductivityInterpolation_ != "unweightedHarmonic"
     && interfaceConductivityInterpolation_ != "distanceWeightedHarmonic"
    )
    {
        FatalErrorInFunction
            << "Unknown interfaceConductivityInterpolation '"
            << interfaceConductivityInterpolation_ << "'. Valid values are "
            << "unweightedHarmonic and distanceWeightedHarmonic."
            << exit(FatalError);
    }
    if
    (
        intracellularAssembly_ != "currentSplit"
     && intracellularAssembly_ != "matchedSubmesh"
    )
    {
        FatalErrorInFunction
            << "Unknown intracellularAssembly '" << intracellularAssembly_
            << "'. Valid values are currentSplit and matchedSubmesh."
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

    if (reportSetup_)
    {
        Info<< "extracellularPotentialDomain: bathZones=" << bathCellZoneNames_
            << " bathConductivityField=" << bathConductivityFieldName_
            << " interfaceConductivityInterpolation="
            << interfaceConductivityInterpolation_
            << " intracellularAssembly=" << intracellularAssembly_
            << " phiENonOrthogonalCorrectors="
            << phiENonOrthogonalCorrectors_;

        if (hasPhiEReferencePoint_)
        {
            Info<< " refPoint=" << phiEReferencePoint_
                << " refValue=" << phiEReferenceValue_;
        }
        else
        {
            Info<< " refPoint=<not required with Dirichlet ground>";
        }

        Info<< endl;
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
    applyGroundPatchValues(phiEPtr_(), dict);

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

    buildHeartScatterMap();
    assembleConductivities();
    buildSigmaTotalSurface();
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

    if (intracellularAssembly_ == "matchedSubmesh")
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
            baseMesh_.time().timeName(),
            baseMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        baseMesh_
    );

    const scalarField& bathSigmaI = bathSigma.primitiveField();
    scalar localMinBath = GREAT;
    scalar localMaxBath = -GREAT;

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
            localMinBath = min(localMinBath, sigma);
            localMaxBath = max(localMaxBath, sigma);
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

    if (reportSetup_)
    {
        reduce(localMinBath, minOp<scalar>());
        reduce(localMaxBath, maxOp<scalar>());
        Info<< "extracellularPotentialDomain: bath sigma min/max = "
            << localMinBath << " / " << localMaxBath << endl;
    }
}

void extracellularPotentialDomain::buildSigmaTotalSurface()
{
    const volTensorField& sigma = sigmaTotalPtr_();
    const volTensorField& sigmaIglobal = sigmaIglobalPtr_();
    const fvMesh& m = baseMesh_;

    sigmaTotalfPtr_.reset
    (
        new surfaceTensorField
        (
            IOobject
            (
                "sigmaTotalf",
                m.time().timeName(),
                m,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            m,
            dimensionedTensor("zero", sigma.dimensions(), tensor::zero)
        )
    );
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

    surfaceTensorField& Sf = sigmaTotalfPtr_();
    surfaceTensorField& Sef = sigmaExtracellularfPtr_();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const tensorField& sigmaI = sigma.primitiveField();
    const tensorField& sigmaIntracellularI = sigmaIglobal.primitiveField();
    const scalarField& weights = m.weights().primitiveField();
    tensorField& SfI = Sf.primitiveFieldRef();
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

        if (interfaceConductivityInterpolation_ == "unweightedHarmonic")
        {
            SfI[faceI] =
                extracellularFaceConductivity::unweightedExtracellularFaceTensor
                (
                    sigmaI[ownCell],
                    sigmaI[neiCell],
                    sigmaIntracellularI[ownCell],
                    sigmaIntracellularI[neiCell]
                );
        }
        else if
        (
            interfaceConductivityInterpolation_
         == "distanceWeightedHarmonic"
        )
        {
            SfI[faceI] =
                extracellularFaceConductivity::
                distanceWeightedExtracellularFaceTensor
                (
                    sigmaI[ownCell],
                    sigmaI[neiCell],
                    sigmaIntracellularI[ownCell],
                    sigmaIntracellularI[neiCell],
                    weights[faceI]
                );
        }
        else
        {
            SfI[faceI] =
                extracellularFaceConductivity::
                distanceWeightedExtracellularFaceTensor
                (
                    sigmaI[ownCell],
                    sigmaI[neiCell],
                    sigmaIntracellularI[ownCell],
                    sigmaIntracellularI[neiCell],
                    weights[faceI]
                );
        }

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

    forAll(Sf.boundaryField(), patchI)
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

            Field<tensor>& Sfp = Sf.boundaryFieldRef()[patchI];
            Field<tensor>& Sefp = Sef.boundaryFieldRef()[patchI];
            const scalarField& patchWeights =
                m.weights().boundaryField()[patchI];

            forAll(Sfp, faceI)
            {
                if
                (
                    interfaceConductivityInterpolation_
                 == "unweightedHarmonic"
                )
                {
                    Sfp[faceI] =
                        extracellularFaceConductivity::
                        unweightedExtracellularFaceTensor
                        (
                            sigmaP[faceI],
                            sigmaN[faceI],
                            sigmaIP[faceI],
                            sigmaIN[faceI]
                        );
                }
                else if
                (
                    interfaceConductivityInterpolation_
                 == "distanceWeightedHarmonic"
                )
                {
                    Sfp[faceI] =
                        extracellularFaceConductivity::
                        distanceWeightedExtracellularFaceTensor
                        (
                            sigmaP[faceI],
                            sigmaN[faceI],
                            sigmaIP[faceI],
                            sigmaIN[faceI],
                            patchWeights[faceI]
                        );
                }
                else
                {
                    Sfp[faceI] = extracellularFaceConductivity::linear
                    (
                        sigmaP[faceI],
                        sigmaN[faceI],
                        patchWeights[faceI]
                    );
                }

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
            // Physical boundaries use the cell-internal value. The sigmaTotal
            // fields use zeroGradient BCs, so patchInternalField is consistent.
            Sf.boundaryFieldRef()[patchI] =
                sigma.boundaryField()[patchI].patchInternalField();
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

    scalarField& heartPhiEI = heartPhiEPtr->primitiveFieldRef();
    const scalarField& globalPhiEI = phiEPtr_().primitiveField();
    forAll(heartCellToBaseCell_, heartCellI)
    {
        heartPhiEI[heartCellI] = globalPhiEI[heartCellToBaseCell_[heartCellI]];
    }
    heartPhiEPtr->correctBoundaryConditions();

    tmp<volScalarField> tHeartRhs = -fvc::laplacian(*GiPtr, *VmPtr);
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

    fvScalarMatrix phiEqn
    (
        fvm::laplacian
        (
            intracellularAssembly_ == "matchedSubmesh"
              ? sigmaExtracellularfPtr_()
              : sigmaTotalfPtr_(),
            phiEPtr_()
        )
     == rhsGlobal
    );

    if (intracellularAssembly_ == "matchedSubmesh")
    {
        fvScalarMatrix heartPhiEqn
        (
            fvm::laplacian(*GiPtr, *heartPhiEPtr)
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
}


void extracellularPotentialDomain::solvePhiE()
{
    scatterHeartVm();

    // Bath/global phiE is a linear elliptic solve with a fixed (scattered-Vm)
    // RHS: it needs the non-orthogonal corrector, not an outer loop. The
    // phiE<->Vm coupling iteration is owned by the advance scheme
    // (bathPdeCouplingMethod), not here.
    correctNonOrthogonalLoop
    (
        phiENonOrthogonalCorrectors_,
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
