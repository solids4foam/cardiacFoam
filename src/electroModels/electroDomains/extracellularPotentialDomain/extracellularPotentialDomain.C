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
#include "compatibilityImbalance.H"
#include "exposedPhiETrace.H"
#include "bidomainSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "OStringStream.H"
#include "dimVoltage.H"
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

// The non-orthogonal corrector count of the PIMPLE dictionary (required)
label resolveNonOrthogonalCorrectors(const fvMesh& baseMesh)
{
    const dictionary& pimpleDict = baseMesh.solutionDict().subDict("PIMPLE");

    const label nCorr = pimpleDict.get<label>("nNonOrthogonalCorrectors");

    if (nCorr < 0)
    {
        FatalIOErrorInFunction(pimpleDict)
            << "nNonOrthogonalCorrectors must be non-negative; found "
            << nCorr << '.' << exit(FatalIOError);
    }

    return nCorr;
}


const volScalarField& requireHeartVm
(
    const myocardiumDomainInterface& heartDomain
)
{
    const volScalarField* VmPtr = heartDomain.VmPtr();

    if (!VmPtr)
    {
        FatalErrorInFunction
            << "Heart domain does not expose Vm."
            << exit(FatalError);
    }

    return *VmPtr;
}


} // End anonymous namespace


extracellularPotentialDomain::extracellularPotentialDomain
(
    const fvMesh& baseMesh,
    myocardiumDomainInterface& heartDomain,
    const dictionary& electroProperties
)
:
    baseMesh_(baseMesh),
    heartDomain_(heartDomain),
    heartVm_(requireHeartVm(heartDomain)),
    conductorPtr_(),
    heartPhiEPtr_(),
    phiEPtr_(),
    VmGlobalPtr_(),
    nNonOrthogonalCorrectors_
    (
        resolveNonOrthogonalCorrectors(baseMesh)
    ),
    manufacturedEvaluationTime_(baseMesh.time().value())
{
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

    // The bidomain intracellular potential is insulated at the heart wall,
    // so the heart side of heart-torso faces carries Ge
    const bool totalWallConductivity = false;

    const dictionary& myocardiumCoeffs = electroProperties.subDict
    (
        electroProperties.get<word>("myocardiumSolver") + "Coeffs"
    );
    const Switch sealedHeartBoundary
    (
        myocardiumCoeffs.get<Switch>("sealedHeartBoundary")
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
            passiveConductor::torsoSettings::read
            (
                baseMesh_, electroProperties.subDict("torso")
            ),
            sealedHeartBoundary,
            heartVm_,
            totalWallConductivity
        )
    );
    const passiveConductor::torsoSettings& torso = conductorPtr_->torso();
    const wordList expectedTypes(torso.patchTypes(baseMesh_, true));

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
            expectedTypes
        )
    );
    {
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

    torso.setGroundValues(phiEPtr_());
    sealedWallPolicy::publish(phiEPtr_(), false);

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

    // The heart phiE, with its final patch types and wall states: the
    // exposed heart faces carry the global trace under bathHeartPhiETrace
    // global (sealed heart only), every heart patch is zeroGradient
    // otherwise
    heartPhiEPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "phiE",
                heartMesh.time().timeName(),
                heartMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            heartMesh,
            dimensionedScalar("phiE", dimVoltage, 0.0),
            "zeroGradient"
        )
    );

    if
    (
        bidomainSolver::readBathHeartPhiETrace(myocardiumCoeffs)
     == bidomainSolver::heartPhiETrace::global
    )
    {
        const label nExposed = useExposedPhiETrace
        (
            heartPhiEPtr_(), *heartFaceMapPtr, baseMesh_
        );

        if (nExposed == 0)
        {
            FatalErrorInFunction
                << "bathHeartPhiETrace global: no exposed heart faces"
                << exit(FatalError);
        }

        Info<< "bathHeartPhiETrace global: " << nExposed
            << " exposed heart faces" << endl;
    }

    sealedWallPolicy::publish(heartPhiEPtr_(), sealedHeartBoundary);

    // The myocardium reads the heart phiE; a verification model writes the
    // manufactured initial global phiE, which the heart phiE then follows
    heartDomain_.bindExternalPhiE
    (
        heartPhiEPtr_(), phiEPtr_(), conductorPtr_->heartCellToBaseCell()
    );
    updateHeartPhiE();
}


void extracellularPotentialDomain::updateHeartPhiE()
{
    const passiveConductor& conductor = conductorPtr_();

    setExposedPhiETrace
    (
        heartPhiEPtr_(), phiEPtr_(), conductor.heartCellToBaseCell(),
        *heartDomain_.subsetFaceMapPtr(),
        *heartDomain_.extracellularConductivityPtr(),
        conductor.sigmaTotal()
    );
}

void extracellularPotentialDomain::scatterHeartVm()
{
    scalarField& VmGlobal = VmGlobalPtr_().primitiveFieldRef();
    const scalarField& heartVm = heartVm_.primitiveField();
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
    const passiveConductor& conductor = conductorPtr_();
    const passiveConductor::torsoSettings& torso = conductor.torso();

    conductor.setSurfaceCurrents(phiEPtr_());

    updateHeartPhiE();

    tmp<volScalarField> trhsGlobal(conductor.heartSource(heartVm_));
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
        conductor.assemble
        (
            phiEPtr_(),
            heartPhiEPtr_(),
            &rhsGlobal,
            storeFlux ? &correction : nullptr
        )
    );
    fvScalarMatrix& phiEqn = tphiEqn.ref();

    if (storeFlux)
    {
        interfaceGradientPtr_() = correction.gradient();
    }

    // Pure Neumann (no ground patch): report the discrete imbalance the
    // reference cell absorbs
    if (torso.pinsReference())
    {
        reportCompatibilityImbalance(phiEqn);
    }

    conductor.applyReference(phiEqn, torso.phiEReferenceValue);

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

    // The myocardium's Vm equation reads the solved potential
    updateHeartPhiE();
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
