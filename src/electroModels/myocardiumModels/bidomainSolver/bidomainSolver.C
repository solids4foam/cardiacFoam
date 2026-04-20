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

defineTypeNameAndDebug(BidomainSolver, 0);
addToRunTimeSelectionTable(myocardiumSolver, BidomainSolver, dictionary);


BidomainSolver::BidomainSolver
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
    phiEReferencePoint_(electroProperties.get<point>("phiERefPoint"))
{
}


label BidomainSolver::referenceCell() const
{
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

tmp<volTensorField> BidomainSolver::initialiseConductivityTensor
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


void BidomainSolver::solveDiffusionExplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt
)
{
    (void)dt;
    const label refCell = referenceCell();

    fvScalarMatrix phiEqn
    (
        fvm::laplacian(GiPlusGe_, phiE_)
     == -fvc::div(Gi_ & fvc::grad(domain.Vm()))
    );
    if (refCell >= 0)
    {
        phiEqn.setReference(refCell, phiEReferenceValue_, true);
    }
    solve(phiEqn);

    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      == fvc::div(Gi_ & fvc::grad(domain.Vm() + phiE_))
       - domain.chi()*domain.Cm()*domain.Iion()
       + domain.sourceField()
    );

    phiI_ = domain.Vm() + phiE_;
}


void BidomainSolver::solveDiffusionImplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt
)
{
    (void)dt;
    const label refCell = referenceCell();

    fvScalarMatrix phiEqn
    (
        fvm::laplacian(GiPlusGe_, phiE_)
     == -fvc::div(Gi_ & fvc::grad(domain.Vm()))
    );
    if (refCell >= 0)
    {
        phiEqn.setReference(refCell, phiEReferenceValue_, true);
    }
    solve(phiEqn);

    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      == fvm::laplacian(Gi_, domain.Vm())
       + fvc::div(Gi_ & fvc::grad(phiE_))
       - domain.chi()*domain.Cm()*domain.Iion()
        + domain.sourceField()
    );

    phiI_ = domain.Vm() + phiE_;
}


void BidomainSolver::solveDiffusionImplicit
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
