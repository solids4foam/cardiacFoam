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
#include "conductivityFieldIO.H"

namespace Foam
{

defineTypeNameAndDebug(monodomainSolver, 0);
addToRunTimeSelectionTable
(
    myocardiumSolver,
    monodomainSolver,
    dictionary
);


monodomainSolver::monodomainSolver
(
    const fvMesh& mesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const dictionary& electroProperties
)
:
    conductivity_
    (
        initialiseConductivity
        (
            mesh,
            supportMesh,
            meshSubsetPtr,
            electroProperties
        )
    )
{}

tmp<volTensorField> monodomainSolver::initialiseConductivity
(
    const fvMesh& mesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const dictionary& electroProperties
) const
{
    return readConductivityField
    (
        mesh,
        supportMesh,
        meshSubsetPtr,
        electroProperties,
        conductivityFieldSpec
        {
            "Conductivity",
            "conductivity",
            "conductivity"
        }
    );
}


void monodomainSolver::solveDiffusionExplicit
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


void monodomainSolver::solveDiffusionImplicit
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

} // End namespace Foam

// ************************************************************************* //
