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

#include "monodomainSolver.H"

#include "IOmanip.H"
#include "myocardiumDomain.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

namespace Foam
{

defineTypeNameAndDebug(MonodomainSolver, 0);
addToRunTimeSelectionTable
(
    myocardiumSolver,
    MonodomainSolver,
    dictionary
);


MonodomainSolver::MonodomainSolver
(
    const fvMesh& mesh,
    const dictionary& electroProperties
)
:
    conductivity_(initialiseConductivity(mesh, electroProperties))
{}

tmp<volTensorField> MonodomainSolver::initialiseConductivity
(
    const fvMesh& mesh,
    const dictionary& electroProperties
) const
{
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                "conductivity",
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
        if (electroProperties.lookupOrDefault<Switch>("reportSetup", false))
        {
            Info << nl
                 << "conductivity not found on disk, using value from "
                 << electroProperties.name()
                 << nl << endl;
        }

        result = dimensionedTensor
        (
            dimensionedSymmTensor
            (
                "conductivity",
                pow3(dimTime) * sqr(dimCurrent)/(dimMass*dimVolume),
                electroProperties
            ) & tensor(I)
        );
    }

    return tresult;
}


void MonodomainSolver::solveDiffusionExplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt
)
{
    (void)dt;

    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      == fvc::laplacian(conductivity_, domain.Vm())
       - domain.chi()*domain.Cm()*domain.Iion()
       + domain.sourceField()
    );
}


void MonodomainSolver::solveDiffusionImplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt
)
{
    (void)dt;

    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      == fvm::laplacian(conductivity_, domain.Vm())
       - domain.chi()*domain.Cm()*domain.Iion()
        + domain.sourceField()
    );
}


void MonodomainSolver::solveDiffusionImplicit
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
