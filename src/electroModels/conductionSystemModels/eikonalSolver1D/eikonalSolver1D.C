#include "eikonalSolver1D.H"
#include "conductionSystemDomain.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(EikonalSolver1D, 0);
addToRunTimeSelectionTable
(
    conductionSystemSolver,
    EikonalSolver1D,
    dictionary
);


EikonalSolver1D::EikonalSolver1D(const dictionary& solverCoeffs)
:
    c0_("c0", solverCoeffs)
{
    Info<< "EikonalSolver1D: conduction velocity c0 = " << c0_.value()
        << " m/s" << endl;
}


void EikonalSolver1D::advance
(
    ConductionSystemDomain& domain,
    scalar t0,
    scalar dt
)
{
    // Eikonal propagation: for each graph edge (i,j) with length L_ij,
    // activation time at j = Tacti + L_ij/c0_.
    // A BFS wavefront sweep is performed until no further updates occur.

    const conductionGraph& G = domain.graph();
    scalarField& Tact = domain.activationTime();

    // Check if at least one node is activated (Tact >= 0)
    if (gMax(Tact) < 0)
    {
        WarningInFunction
            << "No nodes in the conduction system are activated (Tact < 0 for all nodes)." << nl
            << "Eikonal propagation will not start." << endl;
        return;
    }

    const scalar c0 = c0_.value();

    bool updated = true;
    while (updated)
    {
        updated = false;
        forAll(G.edgeNodeA, edgeI)
        {
            label i = G.edgeNodeA[edgeI];
            label j = G.edgeNodeB[edgeI];
            scalar L = G.edgeLengths[edgeI];

            // Forward propagation i -> j
            if (Tact[i] >= 0 && (Tact[j] < 0 || Tact[j] > Tact[i] + L/c0))
            {
                Tact[j] = Tact[i] + L/c0;
                updated = true;
            }
            // Backward propagation j -> i (undirected graph)
            if (Tact[j] >= 0 && (Tact[i] < 0 || Tact[i] > Tact[j] + L/c0))
            {
                Tact[i] = Tact[j] + L/c0;
                updated = true;
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
