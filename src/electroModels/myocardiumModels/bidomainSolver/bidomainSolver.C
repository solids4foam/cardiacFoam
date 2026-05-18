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
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "Switch.H"

namespace Foam
{

defineTypeNameAndDebug(bidomainSolver, 0);
addToRunTimeSelectionTable(myocardiumSolver, bidomainSolver, dictionary);


bidomainSolver::bidomainSolver
(
    const fvMesh& mesh,
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
        "zeroGradient"
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
        word("conductivityIntracellular"),
        electroProperties
    )),
    Ge_(initialiseConductivityTensor
    (
        mesh,
        word("conductivityExtracellular"),
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
    externalPhiECellMapPtr_(nullptr)
{
}


label bidomainSolver::referenceCell() const
{
    if (!hasPhiEReferencePoint_)
    {
        FatalErrorInFunction
            << "bidomainSolver requires phiERefPoint when it solves the "
            << "local extracellular potential. For bidomain-bath cases, "
            << "configure potentialDomain with type "
            << "extracellularPotentialDomain so a global phiE field is bound."
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

    forAll(cellMap, cellI)
    {
        localPhiE[cellI] = globalPhiE[cellMap[cellI]];
    }

    phiE_.correctBoundaryConditions();
}


tmp<volTensorField> bidomainSolver::initialiseConductivityTensor
(
    const fvMesh& mesh,
    const word& fieldName,
    const dictionary& dict
) const
{
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor
            (
                "zero",
                pow3(dimTime) * sqr(dimCurrent)/(dimMass*dimVolume),
                tensor::zero
            )
        )
    );

    volTensorField& result = tresult.ref();

    if (!result.headerOk())
    {
        if (dict.lookupOrDefault<Switch>("reportSetup", false))
        {
            Info << nl
                 << fieldName << " not found on disk, using value from "
                 << dict.name()
                 << nl << endl;
        }

        result = dimensionedTensor
        (
            dimensionedSymmTensor
            (
                fieldName,
                pow3(dimTime) * sqr(dimCurrent)/(dimMass*dimVolume),
                dict
            ) & tensor(I)
        );
    }

    return tresult;
}


void bidomainSolver::solveDiffusionExplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt
)
{
    (void)dt;

    if (externalPhiEBound())
    {
        restrictExternalPhiE();
    }
    else
    {
        const label refCell = referenceCell();

        fvScalarMatrix phiEqn
        (
            fvm::laplacian(GiPlusGe_, phiE_)
         == -fvc::laplacian(Gi_, domain.Vm())
        );
        if (refCell >= 0)
        {
            phiEqn.setReference(refCell, phiEReferenceValue_, true);
        }
        solve(phiEqn);
    }

    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      == fvc::laplacian(Gi_, domain.Vm() + phiE_)
       - domain.chi()*domain.Cm()*domain.Iion()
       + domain.sourceField()
    );

    phiI_ = domain.Vm() + phiE_;
}


void bidomainSolver::solveDiffusionImplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt
)
{
    (void)dt;

    if (externalPhiEBound())
    {
        restrictExternalPhiE();
    }
    else
    {
        const label refCell = referenceCell();

        fvScalarMatrix phiEqn
        (
            fvm::laplacian(GiPlusGe_, phiE_)
         == -fvc::laplacian(Gi_, domain.Vm())
        );
        if (refCell >= 0)
        {
            phiEqn.setReference(refCell, phiEReferenceValue_, true);
        }
        solve(phiEqn);
    }

    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      == fvm::laplacian(Gi_, domain.Vm())
       + fvc::laplacian(Gi_, phiE_)
       - domain.chi()*domain.Cm()*domain.Iion()
        + domain.sourceField()
    );

    phiI_ = domain.Vm() + phiE_;
}


void bidomainSolver::solveDiffusionImplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt,
    pimpleControl& pimple
)
{
    while (pimple.loop())
    {
        solveDiffusionImplicit(domain, dt);
    }
}

} // End namespace Foam

// ************************************************************************* //
