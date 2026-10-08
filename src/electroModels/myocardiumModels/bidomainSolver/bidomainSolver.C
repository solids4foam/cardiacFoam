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
#include "exposedPhiETrace.H"
#include "conormalZeroFluxFvPatchScalarField.H"
#include "fvMeshSubset.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"

namespace Foam
{

defineTypeNameAndDebug(bidomainSolver, 0);
addToRunTimeSelectionTable(myocardiumSolver, bidomainSolver, dictionary);


bidomainSolver::bidomainSolver
(
    const fvMesh& mesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const dictionary& electroProperties
)
:
    phiE_
    (
        IOobject
        (
            "phiE",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiE", dimVoltage, 0.0),
        conormalWallPatchTypes
        (
            mesh,
            electroProperties.get<word>("sealedWallTrace")
        )
    ),
    gradPhiE_
    (
        IOobject
        (
            "grad(" + phiE_.name() + ")",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVoltage/dimLength, vector::zero)
    ),
    phiI_
    (
        IOobject
        (
            "phiI",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
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
    phiEReferenceValue_
    (
        electroProperties.lookupOrDefault<scalar>("phiEReferenceValue", 0.0)
    ),
    phiEReferencePoint_
    (
        electroProperties.found("phiERefPoint")
      ? electroProperties.get<point>("phiERefPoint")
      : point::zero
    ),
    hasPhiEReferencePoint_(electroProperties.found("phiERefPoint")),
    externalPhiEBasePtr_(nullptr),
    externalPhiECellMapPtr_(nullptr),
    meshSubsetPtr_(meshSubsetPtr),
    bathHeartPhiETrace_
    (
        electroProperties.getOrDefault<word>("bathHeartPhiETrace", "zeroGradient")
    )
{
    if
    (
        bathHeartPhiETrace_ != "zeroGradient"
     && bathHeartPhiETrace_ != "global"
    )
    {
        FatalIOErrorInFunction(electroProperties)
            << "bathHeartPhiETrace must be zeroGradient or global"
            << exit(FatalIOError);
    }
    if (bathHeartPhiETrace_ == "global" && !sealedHeartBoundary_)
    {
        FatalIOErrorInFunction(electroProperties)
            << "bathHeartPhiETrace global requires sealedHeartBoundary true"
            << exit(FatalIOError);
    }
    setConormalWallConductivity(phiE_, Ge_.name(), word::null);
}


void bidomainSolver::updateGradPhiE()
{
    gradPhiE_ = fvc::grad(phiE_);
}


label bidomainSolver::referenceCell() const
{
    if (!hasPhiEReferencePoint_)
    {
        FatalErrorInFunction
            << "bidomainSolver requires phiERefPoint when it solves the "
            << "local extracellular potential. For bidomain-bath cases, "
            << "configure bathPotentialDomain so a global phiE field is bound."
            << exit(FatalError);
    }

    const label refCell = phiE_.mesh().findCell(phiEReferencePoint_);

    if (Pstream::parRun())
    {
        const label localOwnsReference = refCell >= 0 ? 1 : 0;
        label ownerCount = localOwnsReference;

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


bool bidomainSolver::externalPhiEBound() const
{
    return externalPhiEBasePtr_ && externalPhiECellMapPtr_;
}


void bidomainSolver::bindExternalPhiE
(
    const volScalarField& phiE,
    const labelUList& heartCellMap
)
{
    externalPhiEBasePtr_ = &phiE;
    externalPhiECellMapPtr_ = &heartCellMap;

    if (bathHeartPhiETrace_ != "global")
    {
        forAll(phiE_.boundaryField(), patchI)
        {
            if
            (
                isA<conormalZeroFluxFvPatchScalarField>
                (
                    phiE_.boundaryField()[patchI]
                )
            )
            {
                FatalErrorInFunction
                    << "sealedWallTrace conormal with a bath requires "
                    << "bathHeartPhiETrace global"
                    << exit(FatalError);
            }
        }
        return;
    }
    if (!meshSubsetPtr_ || !meshSubsetPtr_->hasSubMesh())
    {
        FatalErrorInFunction
            << "bathHeartPhiETrace global requires a myocardium cellZone"
            << exit(FatalError);
    }
    const label nExposed = useExposedPhiETrace
    (
        phiE_, meshSubsetPtr_->faceMap(), phiE.mesh()
    );
    if (nExposed == 0)
    {
        FatalErrorInFunction
            << "bathHeartPhiETrace global: no exposed heart faces"
            << exit(FatalError);
    }
    Info<< "bathHeartPhiETrace global: " << nExposed
        << " exposed heart faces" << endl;
    setExposedPhiETrace
    (
        phiE_, phiE, heartCellMap, meshSubsetPtr_->faceMap(), Ge_,
        phiE.mesh().lookupObject<volTensorField>("sigmaTotal")
    );
}


void bidomainSolver::unbindExternalPhiE()
{
    externalPhiEBasePtr_ = nullptr;
    externalPhiECellMapPtr_ = nullptr;
}


void bidomainSolver::restrictExternalPhiE()
{
    if (!externalPhiEBound())
    {
        FatalErrorInFunction
            << "restrictExternalPhiE called before external phiE was bound."
            << exit(FatalError);
    }

    const scalarField& globalPhiE = externalPhiEBasePtr_->primitiveField();
    const labelUList& cellMap = *externalPhiECellMapPtr_;
    scalarField& localPhiE = phiE_.primitiveFieldRef();

    if (cellMap.size() != localPhiE.size())
    {
        FatalErrorInFunction
            << "External phiE cell map size " << cellMap.size()
            << " does not match bidomain submesh cell count "
            << localPhiE.size() << "."
            << exit(FatalError);
    }

    if
    (
        bathHeartPhiETrace_ == "global"
     && meshSubsetPtr_
     && meshSubsetPtr_->hasSubMesh()
    )
    {
        setExposedPhiETrace
        (
            phiE_, *externalPhiEBasePtr_, cellMap, meshSubsetPtr_->faceMap(),
            Ge_,
            externalPhiEBasePtr_->mesh().lookupObject<volTensorField>("sigmaTotal")
        );
    }
    else
    {
        forAll(cellMap, cellI)
        {
            localPhiE[cellI] = globalPhiE[cellMap[cellI]];
        }
        phiE_.correctBoundaryConditions();
    }
    updateGradPhiE();
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

    if (externalPhiEBound())
    {
        restrictExternalPhiE();
    }
    else
    {
        const label refCell = referenceCell();

        updateGradPhiE();
        fvScalarMatrix phiEqn
        (
            fvm::laplacian
            (
                insulatedFaceConductivity(GiPlusGe_, phiE_, sealedHeartBoundary_),
                phiE_
            )
         == -fvc::laplacian
            (
                insulatedFaceConductivity(Gi_, domain.Vm(), sealedHeartBoundary_),
                domain.Vm()
            )
        );
        if (refCell >= 0)
        {
            phiEqn.setReference(refCell, phiEReferenceValue_, true);
        }
        solve(phiEqn);
    }

    if (const volScalarField* coeff = domain.implicitSourceCoeffPtr())
    {
        solve
        (
            domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
          == fvc::laplacian(tGif(), domain.Vm() + phiE_)
           - domain.chi()*domain.Cm()*domain.Iion()
           + domain.sourceField()
           - (*coeff)*domain.Vm()
        );
    }
    else
    {
        solve
        (
            domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
          == fvc::laplacian(tGif(), domain.Vm() + phiE_)
           - domain.chi()*domain.Cm()*domain.Iion()
           + domain.sourceField()
        );
    }

    phiI_ = domain.Vm() + phiE_;
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
    phiI_ = domain.Vm() + phiE_;
}


void bidomainSolver::solvePhiEImplicitOnce
(
    electroVolumeFieldDomain& domain
)
{
    if (externalPhiEBound())
    {
        restrictExternalPhiE();
    }
    else
    {
        const label refCell = referenceCell();

        updateGradPhiE();
        fvScalarMatrix phiEqn
        (
            fvm::laplacian
            (
                insulatedFaceConductivity(GiPlusGe_, phiE_, sealedHeartBoundary_),
                phiE_
            )
         == -fvc::laplacian
            (
                insulatedFaceConductivity(Gi_, domain.Vm(), sealedHeartBoundary_),
                domain.Vm()
            )
        );
        if (refCell >= 0)
        {
            phiEqn.setReference(refCell, phiEReferenceValue_, true);
        }
        solve(phiEqn);
    }
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

    tmp<volScalarField> tIionExtrap;
    const volScalarField* IionOldPtr = domain.IionOldPtr();
    const volScalarField* IionOldOldPtr = domain.IionOldOldPtr();
    if (IionOldPtr && IionOldOldPtr)
    {
        const scalar deltaT = domain.mesh().time().deltaTValue();
        const scalar deltaT0 = domain.mesh().time().deltaT0Value();
        const scalar r = (deltaT0 > VSMALL) ? (deltaT / deltaT0) : 1.0;
        tIionExtrap = *IionOldPtr + r*(*IionOldPtr - *IionOldOldPtr);
    }
    else
    {
        tIionExtrap = domain.Iion();
    }

    if (const volScalarField* coeff = domain.implicitSourceCoeffPtr())
    {
        solve
        (
            domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
          + fvm::Sp(*coeff, domain.VmRef())
          == fvm::laplacian(tGif(), domain.Vm())
           + fvc::laplacian(tGif(), phiE_)
           - domain.chi()*domain.Cm()*tIionExtrap()
           + domain.sourceField()
        );
    }
    else
    {
        solve
        (
            domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
          == fvm::laplacian(tGif(), domain.Vm())
           + fvc::laplacian(tGif(), phiE_)
           - domain.chi()*domain.Cm()*tIionExtrap()
           + domain.sourceField()
        );
    }
}

} // End namespace Foam

// ************************************************************************* //
