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
    //
    // Finite-volume cable equation on a graph:
    //
    //   dV_i/dt =
    //       (1/(chi*Cm*controlLength_i))
    //         * sum_e sigma_e/L_e * (V_j - V_i)
    //     - Iion_i + Iapp_i/(chi*Cm)
    //
    // where controlLength_i is the 1-D control volume measure assembled from
    // half of each incident edge length. For uniform spacing this reduces to
    // the usual sigma*(V_{i-1}-2V_i+V_{i+1})/dx^2 term.
    const scalar chiCm = domain.chi() * domain.Cm();
    
    // Arrays for the tree solver
    scalarField diag(N, 1.0);       // Main diagonal
    scalarField rhs(N, Zero);       // Right-hand side
    scalarField parentCoeff(N, Zero); // Matrix row child -> column parent
    scalarField childCoeff(N, Zero);  // Matrix row parent -> column child
    scalarField controlLength(N, Zero);

    const labelList& edgeA = domain.edgeStartNodes();
    const labelList& edgeB = domain.edgeEndNodes();
    const scalarField& edgeLength = domain.edgeLengths();
    const scalarField& edgeConductance = domain.edgeConductances();
    
    // Access the topology we built in conductionGraph
    const labelList& parent = domain.graph().parentList;
    const labelList& reverseOrder = domain.graph().reverseOrder;
    const labelList& forwardOrder = domain.graph().orderList;

    forAll(edgeA, edgeI)
    {
        label nodeA = edgeA[edgeI];
        label nodeB = edgeB[edgeI];

        controlLength[nodeA] += 0.5*edgeLength[edgeI];
        controlLength[nodeB] += 0.5*edgeLength[edgeI];
    }

    forAll(controlLength, nodeI)
    {
        if (controlLength[nodeI] <= SMALL)
        {
            FatalErrorInFunction
                << "Node " << nodeI
                << " has zero 1D control length. The cable discretisation "
                << "requires every node to be connected to at least one edge."
                << exit(FatalError);
        }
    }

    // 3a. Assemble finite-volume matrix coefficients
    forAll(edgeA, edgeI)
    {
        label nodeA = edgeA[edgeI];
        label nodeB = edgeB[edgeI];

        const scalar edgeCoeff =
            edgeConductance[edgeI]/edgeLength[edgeI];

        const scalar coeffA =
            dt*edgeCoeff/(chiCm*controlLength[nodeA]);

        const scalar coeffB =
            dt*edgeCoeff/(chiCm*controlLength[nodeB]);

        diag[nodeA] += coeffA;
        diag[nodeB] += coeffB;

        // Determine which node is the child to store the off-diagonal parent link
        if (parent[nodeA] == nodeB)
        {
            parentCoeff[nodeA] = -coeffA;
            childCoeff[nodeA] = -coeffB;
        }
        else if (parent[nodeB] == nodeA)
        {
            parentCoeff[nodeB] = -coeffB;
            childCoeff[nodeB] = -coeffA;
        }
        else
        {
            FatalErrorInFunction
                << "Edge " << edgeI << " connecting nodes " << nodeA
                << " and " << nodeB
                << " is inconsistent with the rooted tree topology."
                << exit(FatalError);
        }
    }

    // 3b. Assemble RHS. The ionic model computes Iion during the reaction
    // sub-step, but Vm is not advanced inside the ODE solver for monodomain
    // runs, so Iion must still enter the membrane balance here.
    for (label i = 0; i < N; i++)
    {
        rhs[i] =
            Vm[i]
          - dt*Iion[i]
          + dt*appliedCurrentBuffer_[i]/chiCm;
    }

    // 3c. Forward Sweep (Bottom-up: Leaves to Root)
    // Iterate through reverseOrder, skipping the last element (which is the root)
    for (label i = 0; i < N - 1; i++)
    {
        label c = reverseOrder[i]; // Child
        label p = parent[c];       // Parent

        scalar factor = childCoeff[c] / diag[c];
        diag[p] -= factor * parentCoeff[c];
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
