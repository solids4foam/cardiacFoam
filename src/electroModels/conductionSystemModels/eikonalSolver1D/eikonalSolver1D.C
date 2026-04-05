#include "eikonalSolver1D.H"
#include "graphConductionSystemDomain.H"
#include "conductionGraph.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(EikonalSolver1D, 0);
addToRunTimeSelectionTable
(
    GraphConductionSystemSolver,
    EikonalSolver1D,
    dictionary
);


EikonalSolver1D::EikonalSolver1D(const dictionary& solverCoeffs)
:
    c0_(solverCoeffs.get<scalar>("c0"))
{
    Info<< "EikonalSolver1D: conduction velocity c0 = " << c0_
        << " m/s" << endl;
}


void EikonalSolver1D::advance
(
    GraphConductionSystemDomain& domain,
    scalar t0,
    scalar dt
)
{
    // Eikonal propagation: for each graph edge (i,j) with length L_ij,
    // activation time at j = Tacti + L_ij/c0_.
    // A BFS wavefront sweep is performed until no further updates occur.

    const conductionGraph& G = domain.graph();
    scalarField& Tact = domain.activationTime();

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
            if (Tact[i] >= 0 && (Tact[j] < 0 || Tact[j] > Tact[i] + L/c0_))
            {
                Tact[j] = Tact[i] + L/c0_;
                updated = true;
            }
            // Backward propagation j -> i (undirected graph)
            if (Tact[j] >= 0 && (Tact[i] < 0 || Tact[i] > Tact[j] + L/c0_))
            {
                Tact[i] = Tact[j] + L/c0_;
                updated = true;
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
