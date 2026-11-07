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

#include "bidomainSolver.H"

#include "IOmanip.H"
#include "PstreamReduceOps.H"
#include "myocardiumDomain.H"
#include "insulatedFaceConductivity.H"
#include "sealedWallPolicy.H"
#include "compatibilityImbalance.H"
#include "conormalStabilisation.H"
#include "fvMeshSubset.H"
#include "torsoCoupling.H"
#include "heartPhiETrace.H"
#include "passiveConductor.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"

namespace Foam
{

defineTypeNameAndDebug(bidomainSolver, 0);
addToRunTimeSelectionTable(myocardiumSolver, bidomainSolver, dictionary);


namespace
{

// With torsoCoupling bath the bath domain owns phiE and binds it to the
// heart; the heart reads no phiE of its own
bool solvedWithBath(const dictionary& electroProperties)
{
    return
        torsoCoupling::readFromCoeffs(electroProperties)
     == torsoCoupling::type::bath;
}


// The reference of the heart's own phiE solve: phiERefPoint and the value
// pinned at its cell, both required unless the bath owns phiE
autoPtr<Tuple2<point, scalar>> readPhiEReference
(
    const dictionary& electroProperties
)
{
    if
    (
        electroProperties.found("phiEReferenceValue")
     && !electroProperties.found("phiERefPoint")
    )
    {
        FatalIOErrorInFunction(electroProperties)
            << "phiEReferenceValue is set but phiERefPoint is not; the value "
            << "is applied at the cell containing phiERefPoint"
            << exit(FatalIOError);
    }

    if (solvedWithBath(electroProperties))
    {
        return nullptr;
    }

    return autoPtr<Tuple2<point, scalar>>::New
    (
        electroProperties.get<point>("phiERefPoint"),
        electroProperties.get<scalar>("phiEReferenceValue")
    );
}

} // End anonymous namespace


bidomainSolver::bidomainSolver
(
    const fvMesh& mesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const dictionary& electroProperties
)
:
    ownPhiEPtr_(),
    phiEPtr_(nullptr),
    phiI_
    (
        IOobject
        (
            "phiI",
            mesh.time().timeName(),
            mesh,
            myocardiumFieldReadOption(meshSubsetPtr),
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiI", dimVoltage, 0.0),
        "zeroGradient"
    ),
    Gi_(initialiseConductivityTensor
    (
        mesh,
        supportMesh,
        meshSubsetPtr,
        conductivityFieldSpec
        {
            "ConductivityIntracellular",
            "conductivityIntracellular",
            "conductivityIntracellular"
        },
        electroProperties
    )),
    Ge_(initialiseConductivityTensor
    (
        mesh,
        supportMesh,
        meshSubsetPtr,
        conductivityFieldSpec
        {
            "ConductivityExtracellular",
            "conductivityExtracellular",
            "conductivityExtracellular"
        },
        electroProperties
    )),
    GiPlusGe_
    (
        IOobject
        (
            "GiPlusGe",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Gi_ + Ge_
    ),
    sealedHeartBoundary_
    (
        electroProperties.get<Switch>("sealedHeartBoundary")
    ),
    phiEReferencePtr_(readPhiEReference(electroProperties)),
    verifierPtr_(nullptr)
{
    // Validates bathHeartPhiETrace: required with a torso, rejected without
    // one; the bath domain reads it for its heart phiE
    (void)heartPhiETrace::read(electroProperties);

    // With a bath the bath domain owns, solves and writes the heart phiE,
    // and binds it here (bindExternalPhiE)
    if (!solvedWithBath(electroProperties))
    {
        ownPhiEPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "phiE",
                    mesh.time().timeName(),
                    mesh,
                    myocardiumFieldReadOption(meshSubsetPtr),
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("phiE", dimVoltage, 0.0),
                "zeroGradient"
            )
        );
        phiEPtr_ = ownPhiEPtr_.get();
        sealedWallPolicy::publish(ownPhiEPtr_(), sealedHeartBoundary_);
    }

    if (meshSubsetPtr && meshSubsetPtr->hasSubMesh())
    {
        const labelUList& cellMap = meshSubsetPtr->cellMap();

        if (ownPhiEPtr_)
        {
            readMappedCellField(ownPhiEPtr_(), supportMesh, cellMap);
        }
        readMappedCellField(phiI_, supportMesh, cellMap);
    }
}


label bidomainSolver::referenceCell() const
{
    if (!phiEReferencePtr_)
    {
        FatalErrorInFunction
            << "bidomainSolver requires phiERefPoint when it solves the "
            << "local extracellular potential. For bidomain-bath cases, "
            << "set torso and torsoCoupling bath so a global phiE field is "
            << "bound."
            << exit(FatalError);
    }

    return passiveConductor::locateReferenceCell
    (
        ownPhiEPtr_().mesh(), phiEReferencePtr_->first()
    );
}


const volScalarField& bidomainSolver::phiE() const
{
    if (!phiEPtr_)
    {
        FatalErrorInFunction
            << "The heart phiE of the bath has not been bound: the bath "
            << "domain binds it when it is constructed."
            << exit(FatalError);
    }

    return *phiEPtr_;
}


void bidomainSolver::bindExternalPhiE(const volScalarField& heartPhiE)
{
    if (ownPhiEPtr_)
    {
        FatalErrorInFunction
            << "This bidomainSolver owns and solves its phiE (no "
            << "torsoCoupling bath); it cannot read another one."
            << exit(FatalError);
    }

    if (&heartPhiE.mesh() != &phiI_.mesh())
    {
        FatalErrorInFunction
            << "The bound phiE lives on mesh '" << heartPhiE.mesh().name()
            << "', not on the myocardium mesh '" << phiI_.mesh().name()
            << "'." << exit(FatalError);
    }

    phiEPtr_ = &heartPhiE;
}


void bidomainSolver::unbindExternalPhiE()
{
    phiEPtr_ = ownPhiEPtr_.get();
}



tmp<volTensorField> bidomainSolver::initialiseConductivityTensor
(
    const fvMesh& mesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const conductivityFieldSpec& spec,
    const dictionary& dict
) const
{
    return readConductivityField
    (
        mesh,
        supportMesh,
        meshSubsetPtr,
        dict,
        spec
    );
}


void bidomainSolver::solveDiffusionExplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt
)
{
    (void)dt;

    const tmp<surfaceTensorField> tGif
    (
        insulatedFaceConductivity(Gi_, domain.Vm(), sealedHeartBoundary_)
    );

    solvePhiEImplicitOnce(domain);

    // The junction coefficient is cell-local, so the matrix stays diagonal
    // and the explicit solve keeps its cost.
    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      + fvm::Sp(domain.implicitSourceCoeff(), domain.VmRef())
     == fvc::laplacian(tGif(), domain.Vm())
      + fvc::laplacian(tGif(), phiE())
      - domain.chi()*domain.Cm()*domain.Iion()
      + domain.sourceField()
    );

    phiI_ = domain.Vm() + phiE();
}


void bidomainSolver::solveDiffusionImplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt
)
{
    (void)dt;

    // Heart-only bidomain: phiE and Vm are re-solved together on every outer
    // PIMPLE corrector (myocardiumSolver::solveDiffusionImplicit). The bath
    // bidomain runs its outer PIMPLE passes in
    // staggeredElectrophysicsAdvanceScheme through solveDiffusionStepOnce.
    solvePhiEImplicitOnce(domain);
    solveVmImplicitOnce(domain);
    phiI_ = domain.Vm() + phiE();
}


void bidomainSolver::solvePhiEImplicitOnce
(
    electroVolumeFieldDomain& domain
)
{
    if (!ownPhiEPtr_)
    {
        // The bath domain solves the global phiE and updates the bound
        // heart phiE; the Vm equation reads it as it is
        (void)phiE();
        return;
    }

    volScalarField& phiE = ownPhiEPtr_();
    const label refCell = referenceCell();

    volScalarField rhs
    (
        -fvc::laplacian
        (
            insulatedFaceConductivity(Gi_, domain.Vm(), sealedHeartBoundary_),
            domain.Vm()
        )
    );
    if (verifierPtr_)
    {
        verifierPtr_->addManufacturedPhiESource(rhs, phiE.time().value());
    }
    const tmp<surfaceTensorField> tGf
    (
        insulatedFaceConductivity(GiPlusGe_, phiE, sealedHeartBoundary_)
    );
    fvScalarMatrix phiEqn
    (
        fvm::laplacian(tGf(), phiE)
      + conormalStabilisation(tGf(), phiE)
     == rhs
    );
    // Pure Neumann (no patch fixes the value): report the discrete
    // imbalance the reference cell absorbs, then pin that cell
    if (phiE.needReference())
    {
        reportCompatibilityImbalance(phiEqn);
        phiEqn.setReference(refCell, phiEReferencePtr_->second());
    }
    solve(phiEqn);
}


void bidomainSolver::solveVmImplicitOnce
(
    electroVolumeFieldDomain& domain
)
{
    const tmp<surfaceTensorField> tGif
    (
        insulatedFaceConductivity(Gi_, domain.Vm(), sealedHeartBoundary_)
    );

    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      + fvm::Sp(domain.implicitSourceCoeff(), domain.VmRef())
     == fvm::laplacian(tGif(), domain.Vm())
      + fvc::laplacian(tGif(), phiE())
      - domain.chi()*domain.Cm()*domain.ionicCurrentForImplicitDiffusion()
      + domain.sourceField()
    );
}

} // End namespace Foam

// ************************************************************************* //
