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


tensor harmonicFaceTensor(const tensor& a, const tensor& b)
{
    tensor result(tensor::zero);

    for (direction d = 0; d < tensor::nComponents; ++d)
    {
        const scalar denom = a[d] + b[d];
        result[d] =
            (mag(denom) > SMALL && a[d]*b[d] > 0.0)
          ? 2.0*a[d]*b[d]/denom
          : 0.5*denom;
    }

    return result;
}


tensor extracellularFaceTensor
(
    const tensor& sP,
    const tensor& sN,
    const tensor& sigmaIP,
    const tensor& sigmaIN
)
{
    const bool pIsHeart = magSqr(sigmaIP) > SMALL;
    const bool nIsHeart = magSqr(sigmaIN) > SMALL;

    if (pIsHeart != nIsHeart)
    {
        const tensor heartSigmaE = pIsHeart ? sP - sigmaIP : sN - sigmaIN;
        const tensor bathSigma = pIsHeart ? sN : sP;

        return harmonicFaceTensor(heartSigmaE, bathSigma);
    }

    return harmonicFaceTensor(sP, sN);
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
    heartCellZoneName_(dict.lookupOrDefault<word>("heartCellZone", "myocardium")),
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
    reportSetup_(dict.lookupOrDefault<Switch>("reportSetup", false))
{
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
        Info<< "extracellularPotentialDomain: heartZone=" << heartCellZoneName_
            << " bathZones=" << bathCellZoneNames_
            << " bathConductivityField=" << bathConductivityFieldName_;

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
    buildHarmonicSigmaTotalSurface();
    heartDomain_.bindExternalPhiE(phiEPtr_(), heartCellToBaseCell_);
}

void extracellularPotentialDomain::buildHeartScatterMap()
{
    const labelUList* heartCellMapPtr = heartDomain_.subsetCellMapPtr();

    if (!heartCellMapPtr)
    {
        FatalErrorInFunction
            << "Unified phiE requires the myocardium domain to expose a "
            << "subset cell map to the base mesh."
            << exit(FatalError);
    }

    heartCellToBaseCell_ = labelList(*heartCellMapPtr);
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

    if (heartCellToBaseCell_.size() != Gi.size())
    {
        FatalErrorInFunction
            << "Heart/base cell map size " << heartCellToBaseCell_.size()
            << " does not match conductivity cell count " << Gi.size() << "."
            << exit(FatalError);
    }

    forAll(heartCellToBaseCell_, heartCellI)
    {
        const label baseCellI = heartCellToBaseCell_[heartCellI];
        sigmaTotalI[baseCellI] = Gi[heartCellI] + Ge[heartCellI];
        sigmaIglobalI[baseCellI] = Gi[heartCellI];
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
            sigmaTotalI[cellI] = sigma*tensor::I;
            localMinBath = min(localMinBath, sigma);
            localMaxBath = max(localMaxBath, sigma);
        }
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

void extracellularPotentialDomain::buildHarmonicSigmaTotalSurface()
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

    surfaceTensorField& Sf = sigmaTotalfPtr_();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const tensorField& sigmaI = sigma.primitiveField();
    const tensorField& sigmaIntracellularI = sigmaIglobal.primitiveField();
    tensorField& SfI = Sf.primitiveFieldRef();

    forAll(own, faceI)
    {
        const label ownCell = own[faceI];
        const label neiCell = nei[faceI];

        SfI[faceI] =
            extracellularFaceTensor
            (
                sigmaI[ownCell],
                sigmaI[neiCell],
                sigmaIntracellularI[ownCell],
                sigmaIntracellularI[neiCell]
            );
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

            forAll(Sfp, faceI)
            {
                Sfp[faceI] =
                    extracellularFaceTensor
                    (
                        sigmaP[faceI],
                        sigmaN[faceI],
                        sigmaIP[faceI],
                        sigmaIN[faceI]
                    );
            }
        }
        else
        {
            // Physical boundaries use the cell-internal value. The sigmaTotal
            // fields use zeroGradient BCs, so patchInternalField is consistent.
            Sf.boundaryFieldRef()[patchI] =
                sigma.boundaryField()[patchI].patchInternalField();
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

    if (Pstream::parRun())
    {
        const label localOwns = refCell >= 0 ? 1 : 0;
        label ownerCount = localOwns;
        reduce(ownerCount, sumOp<label>());

        if (ownerCount != 1)
        {
            FatalErrorInFunction
                << "phiERefPoint " << phiEReferencePoint_
                << " must be owned by exactly one processor in parallel, but "
                << ownerCount << " processors reported a containing cell."
                << exit(FatalError);
        }
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

void extracellularPotentialDomain::solvePhiE()
{
    scatterHeartVm();

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
        fvm::laplacian(sigmaTotalfPtr_(), phiEPtr_())
     == rhsGlobal
    );

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

void extracellularPotentialDomain::prepareTimeStep(scalar t0, scalar dt)
{
    (void)t0; (void)dt;
    scatterHeartVm();
}

void extracellularPotentialDomain::advance(scalar t0, scalar dt)
{
    (void)t0; (void)dt;
    scatterHeartVm();
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
