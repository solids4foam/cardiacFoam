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
#include "conductionSystemDomain.H"
#include "ionicModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(Monodomain1DSolver, 0);
addToRunTimeSelectionTable
(
    conductionSystemSolver,
    Monodomain1DSolver,
    dictionary
);

void Monodomain1DSolver::advance
(
    ConductionSystemDomain& domain,
    scalar t0,
    scalar dt
)
{
    scalarField& Vm = domain.membranePotential();
    scalarField& Iion = domain.ionicCurrent();
    const label N = Vm.size();

    // ---- Step 1: Reaction (First Half-Step) ----
    // Integrate ODEs for dt/2
    domain.ionicModelRef().solveODE(t0, dt / 2.0, Vm, Iion);

    // ---- Step 2: Applied current ----
    appliedCurrentBuffer_.setSize(N);
    domain.assembleAppliedCurrent(t0, appliedCurrentBuffer_);

    // ---- Step 3: Diffusion (Full-Step PDE) using O(N) Hines Algorithm ----
    const scalar dtOverChiCm = dt / (domain.chi() * domain.Cm());
    
    // Arrays for the tree solver
    scalarField diag(N, 1.0);       // Main diagonal
    scalarField rhs(N, Zero);       // Right-hand side
    scalarField parentCoeff(N, Zero); // Off-diagonal coefficient linking to parent

    const labelList& edgeA = domain.edgeStartNodes();
    const labelList& edgeB = domain.edgeEndNodes();
    const scalarField& edgeLength = domain.edgeLengths();
    const scalarField& edgeConductance = domain.edgeConductances();
    
    // Access the topology we built in conductionGraph
    const labelList& parent = domain.graph().parentList;
    const labelList& reverseOrder = domain.graph().reverseOrder;
    const labelList& forwardOrder = domain.graph().orderList;

    // 3a. Assemble matrix coefficients
    forAll(edgeA, edgeI)
    {
        label nodeA = edgeA[edgeI];
        label nodeB = edgeB[edgeI];
        
        scalar coeff = dtOverChiCm * edgeConductance[edgeI] / sqr(edgeLength[edgeI]);

        diag[nodeA] += coeff;
        diag[nodeB] += coeff;

        // Determine which node is the child to store the off-diagonal parent link
        if (parent[nodeA] == nodeB)
        {
            parentCoeff[nodeA] = -coeff;
        }
        else if (parent[nodeB] == nodeA)
        {
            parentCoeff[nodeB] = -coeff;
        }
    }

    // 3b. Assemble RHS (Notice Iion is REMOVED here because ODE handles it)
    for (label i = 0; i < N; i++)
    {
        rhs[i] = Vm[i] + dtOverChiCm * appliedCurrentBuffer_[i];
    }

    // 3c. Forward Sweep (Bottom-up: Leaves to Root)
    // Iterate through reverseOrder, skipping the last element (which is the root)
    for (label i = 0; i < N - 1; i++)
    {
        label c = reverseOrder[i]; // Child
        label p = parent[c];       // Parent

        scalar factor = parentCoeff[c] / diag[c];
        diag[p] -= factor * parentCoeff[c]; // Because symmetric, child-to-parent is same as parent-to-child
        rhs[p]  -= factor * rhs[c];
    }

    // 3d. Backward Sweep (Top-down: Root to Leaves)
    // Solve the root node first
    label root = forwardOrder[0];
    Vm[root] = rhs[root] / diag[root];

    // Propagate the exact solution down to the leaves
    for (label i = 1; i < N; i++)
    {
        label c = forwardOrder[i];
        label p = parent[c];

        Vm[c] = (rhs[c] - parentCoeff[c] * Vm[p]) / diag[c];
    }

    // ---- Step 4: Reaction (Second Half-Step) ----
    // Integrate ODEs for the remaining dt/2 using the newly diffused Vm
    domain.ionicModelRef().solveODE(t0 + (dt / 2.0), dt / 2.0, Vm, Iion);

    domain.reportAdvanceDiagnostics(t0, dt);
}

} // End namespace Foam 
