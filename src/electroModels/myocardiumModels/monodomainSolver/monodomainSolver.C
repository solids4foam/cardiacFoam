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
#include "electroVolumeFieldDomain.H"
#include "dimVoltage.H"
#include "addToRunTimeSelectionTable.H"
#include "conductivityFieldIO.H"
#include "insulatedFaceConductivity.H"

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
    const dictionary& electroProperties,
    const bool externalPhiE
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
    ),
    sealedHeartBoundary_
    (
        electroProperties.get<Switch>("sealedHeartBoundary")
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

    const tmp<surfaceTensorField> tGf
    (
        insulatedFaceConductivity
        (
            conductivity_, domain.Vm(), sealedHeartBoundary_
        )
    );

    // The junction coefficient is cell-local, so the matrix stays diagonal
    // and the explicit solve keeps its cost.
    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      + fvm::Sp(domain.implicitSourceCoeff(), domain.VmRef())
     == fvc::laplacian(tGf(), domain.Vm())
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

    const tmp<surfaceTensorField> tGf
    (
        insulatedFaceConductivity
        (
            conductivity_, domain.Vm(), sealedHeartBoundary_
        )
    );

    solve
    (
        domain.chi()*domain.Cm()*fvm::ddt(domain.VmRef())
      + fvm::Sp(domain.implicitSourceCoeff(), domain.VmRef())
     == fvm::laplacian(tGf(), domain.Vm())
      - domain.chi()*domain.Cm()*domain.ionicCurrentForImplicitDiffusion()
      + domain.sourceField()
    );
}

} // End namespace Foam

// ************************************************************************* //
