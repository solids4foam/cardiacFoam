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
#include "electroVerificationModel.H"
#include "sealedWallPolicy.H"
#include "compatibilityProjection.H"
#include "exposedPhiETrace.H"
#include "addToRunTimeSelectionTable.H"
#include "OStringStream.H"
#include "dimVoltage.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"
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


wordList potentialPatchTypes
(
    const fvMesh& mesh,
    const passiveConductor::torsoSettings& torso
)
{
    wordList patchTypes
    (
        mesh.boundary().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(torso.groundPatches, patchI)
    {
        const label patchID =
            mesh.boundaryMesh().findPatchID(torso.groundPatches[patchI]);
        patchTypes[patchID] = fixedValueFvPatchScalarField::typeName;
    }

    forAll(torso.surfaceCurrentPatches, patchI)
    {
        patchTypes
        [
            mesh.boundaryMesh().findPatchID(torso.surfaceCurrentPatches[patchI])
        ] = fixedGradientFvPatchScalarField::typeName;
    }

    return patchTypes;
}


void applyGroundPatchValues
(
    volScalarField& phiE,
    const passiveConductor::torsoSettings& torso
)
{
    forAll(torso.groundPatches, patchI)
    {
        const label patchId =
            phiE.mesh().boundaryMesh().findPatchID(torso.groundPatches[patchI]);

        phiE.boundaryFieldRef()[patchId] = torso.groundValues[patchI];
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
    torso_(passiveConductor::torsoSettings::read(baseMesh, dict)),
    phiEPtr_(),
    VmGlobalPtr_(),
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
            potentialPatchTypes(baseMesh_, torso_)
        )
    );
    {
        const wordList expectedTypes(potentialPatchTypes(baseMesh_, torso_));
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

    applyGroundPatchValues(phiEPtr_(), torso_);
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

    const labelUList* heartCellMapPtr = heartDomain_.subsetCellMapPtr();
    const labelUList* heartFaceMapPtr = heartDomain_.subsetFaceMapPtr();
    const fvMesh& heartMesh =
        static_cast<const electroStateProvider&>(heartDomain_).mesh();

    if (!heartCellMapPtr)
    {
        FatalErrorInFunction
            << "Unified phiE requires the myocardium domain to expose a "
            << "subset cell map to the base mesh."
            << exit(FatalError);
    }

    if (!heartFaceMapPtr)
    {
        FatalErrorInFunction
            << "Matched intracellular assembly requires the heart/base "
            << "face map."
            << exit(FatalError);
    }

    const volTensorField* GiPtr = heartDomain_.intracellularConductivityPtr();
    const volTensorField* GePtr = heartDomain_.extracellularConductivityPtr();

    if (!GiPtr || !GePtr)
    {
        FatalErrorInFunction
            << "extracellularPotentialDomain requires split intracellular and "
            << "extracellular conductivities from a bidomain myocardium."
            << exit(FatalError);
    }

    const volScalarField* heartVmPtr = heartDomain_.VmPtr();

    if (!heartVmPtr)
    {
        FatalErrorInFunction
            << "Heart domain does not expose Vm."
            << exit(FatalError);
    }

    const autoPtr<volScalarField> bathSigmaPtr
    (
        passiveConductor::readTorsoConductivity
        (
            baseMesh_, torso_.bathConductivityField, true
        )
    );

    conductorPtr_.reset
    (
        new passiveConductor
        (
            baseMesh_,
            heartMesh,
            *heartCellMapPtr,
            *heartFaceMapPtr,
            *GiPtr,
            *GePtr,
            bathSigmaPtr(),
            torso_.bathCellZones,
            torso_.interfaceConductivityInterpolation,
            sealedHeartBoundary_,
            *heartVmPtr,
            false
        )
    );
    heartDomain_.bindExternalPhiE
    (
        phiEPtr_(), conductorPtr_->heartCellToBaseCell()
    );
}

label extracellularPotentialDomain::referenceCell() const
{
    torso_.requirePhiERefPoint("extracellularPotentialDomain");

    return conductorPtr_->referenceCell(torso_.phiERefPoint);
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

    const labelList& heartCellToBaseCell =
        conductorPtr_->heartCellToBaseCell();

    forAll(heartCellToBaseCell, heartCellI)
    {
        VmGlobal[heartCellToBaseCell[heartCellI]] = heartVm[heartCellI];
    }

    VmGlobalPtr_().correctBoundaryConditions();
}

const electroVerificationModel* extracellularPotentialDomain::verifier() const
{
    return heartDomain_.verificationModel();
}


void extracellularPotentialDomain::storeInterfaceFlux
(
    const surfaceScalarField& F
)
{
    surfaceScalarField& flux = interfaceFluxPtr_();
    flux = dimensionedScalar(flux.dimensions(), Zero);

    const volTensorField& sigmaI = conductorPtr_->sigmaIntracellular();
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
}


void extracellularPotentialDomain::solvePhiEOnce()
{
    forAll(torso_.surfaceCurrentPatches, patchI)
    {
        const label patchId =
            baseMesh_.boundaryMesh().findPatchID(torso_.surfaceCurrentPatches[patchI]);

        if (patchId < 0)
        {
            FatalErrorInFunction
                << "Cannot find extracellularPotentialDomain surface-current patch '"
                << torso_.surfaceCurrentPatches[patchI] << "' on mesh '"
                << baseMesh_.name() << "'."
                << exit(FatalError);
        }

        fixedGradientFvPatchScalarField& phiEPatch =
            refCast<fixedGradientFvPatchScalarField>
            (
                phiEPtr_().boundaryFieldRef()[patchId]
            );

        tmp<tensorField> tSigmaPatch =
            conductorPtr_->sigmaTotal()
                .boundaryField()[patchId].patchInternalField();
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
                    << torso_.surfaceCurrentPatches[patchI]
                    << "' found non-positive normal conductivity " << sigmaN
                    << " on face " << faceI << "."
                    << exit(FatalError);
            }

            gradient[faceI] = torso_.surfaceCurrentValues[patchI]/sigmaN;
        }

        phiEPatch.evaluate();
    }

    const volScalarField* VmPtr = heartDomain_.VmPtr();

    if (!VmPtr)
    {
        FatalErrorInFunction
            << "extracellularPotentialDomain requires the heart Vm field."
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
        *heartPhiEPtr, phiEPtr_(), conductorPtr_->heartCellToBaseCell(),
        *heartDomain_.subsetFaceMapPtr(),
        *heartDomain_.extracellularConductivityPtr(),
        conductorPtr_->sigmaTotal()
    );

    tmp<volScalarField> trhsGlobal(conductorPtr_->heartSource(*VmPtr));
    volScalarField& rhsGlobal = trhsGlobal.ref();

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
    passiveConductor::faceCorrection correction;

    tmp<fvScalarMatrix> tphiEqn
    (
        conductorPtr_->assemble
        (
            phiEPtr_(),
            *heartPhiEPtr,
            &rhsGlobal,
            storeFlux ? &correction : nullptr
        )
    );
    fvScalarMatrix& phiEqn = tphiEqn.ref();

    if (storeFlux)
    {
        interfaceGradientPtr_() = correction.gradient();
    }

    if
    (
        !hasDirichletPatch_
     && verifier()
     && verifier()->requiresPhiECompatibilityProjection()
    )
    {
        projectCompatibility(phiEqn);
    }

    if (!hasDirichletPatch_)
    {
        conductorPtr_->applyReference
        (
            phiEqn, phiEPtr_(), referenceCell(), torso_.phiEReferenceValue
        );
    }

    solve(phiEqn);

    if (storeFlux)
    {
        tmp<surfaceScalarField> tF(phiEqn.flux());
        if (correction.crossFlux)
        {
            tF.ref() += correction.crossFlux();
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
    if (conductorPtr_.valid()) conductorPtr_->sigmaTotal().write();
    if (VmGlobalPtr_.valid()) VmGlobalPtr_().write();
}

} // End namespace Foam

// ************************************************************************* //
