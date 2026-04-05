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

#include "monodomain1DSolver.H"
#include "graphConductionSystemDomain.H"
#include "ionicModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(Monodomain1DSolver, 0);
addToRunTimeSelectionTable
(
    GraphConductionSystemSolver,
    Monodomain1DSolver,
    dictionary
);


void Monodomain1DSolver::advance
(
    GraphConductionSystemDomain& domain,
    scalar t0,
    scalar dt
)
{
    scalarField& Vm = domain.membranePotential();
    scalarField& Iion = domain.ionicCurrent();
    const label N = Vm.size();

    // ---- Step 1: Reaction — solve ionic model ODE ----
    domain.ionicModelRef().solveODE(t0, dt, Vm, Iion);

    // ---- Step 2: Applied current (root stimulus + PVJ coupling) ----
    appliedCurrentBuffer_.setSize(N);
    domain.assembleAppliedCurrent(t0, appliedCurrentBuffer_);

    // ---- Step 3: Diffusion — implicit cable equation ----
    //
    //  Point-wise cable equation (backward Euler):
    //
    //    (Vm_new - Vm_old) / dt
    //        = 1/(chi*Cm) * sum_edges[ sigma/L^2 * (Vm_new_j - Vm_new_i) ]
    //          - Iion
    //          + Iapp / (chi*Cm)
    //
    //  Rearranged as  A * Vm_new = b :
    //
    //    A[i][i] = 1 + dt/(chi*Cm) * sum_edges_at_i( sigma_e / L_e^2 )
    //    A[i][j] =   - dt/(chi*Cm) * sigma_e / L_e^2
    //    b[i]    = Vm_old[i] - dt*Iion[i] + dt/(chi*Cm)*Iapp[i]
    //
    //  This is the same implicit discretisation that
    //  fvm::ddt + fvm::laplacian produces on the 3-D myocardium mesh.

    scalarSquareMatrix A(N, Zero);
    scalarField rhs(N, Zero);

    const scalar dtOverChiCm = dt / (domain.chi() * domain.Cm());

    // Diagonal: identity from backward-Euler ddt
    for (label i = 0; i < N; i++)
    {
        A(i, i) = 1.0;
    }

    // Implicit Laplacian contributions from graph edges
    const labelList& edgeA = domain.edgeStartNodes();
    const labelList& edgeB = domain.edgeEndNodes();
    const scalarField& edgeLength = domain.edgeLengths();
    const scalarField& edgeConductance = domain.edgeConductances();

    forAll(edgeA, edgeI)
    {
        const label iN = edgeA[edgeI];
        const label jN = edgeB[edgeI];
        const scalar coeff =
            dtOverChiCm
          * edgeConductance[edgeI]
          / sqr(edgeLength[edgeI]);

        A(iN, iN) += coeff;
        A(jN, jN) += coeff;
        A(iN, jN) -= coeff;
        A(jN, iN) -= coeff;
    }

    // Right-hand side
    forAll(Vm, i)
    {
        rhs[i] = Vm[i]
               - dt * Iion[i]
               + dtOverChiCm * appliedCurrentBuffer_[i];
    }

    // ---- Step 4: Solve  A * x = b  (LU decomposition) ----
    LUsolve(A, rhs);

    // ---- Step 5: Update Vm (rhs now contains the solution) ----
    Vm = rhs;

    domain.reportAdvanceDiagnostics(t0, dt);
}

} // End namespace Foam

// ************************************************************************* //
