#include "eikonalSolver1D.H"
#include "conductionSystemDomain.H"
#include "addToRunTimeSelectionTable.H"
#include "DynamicList.H"
#include <queue>
#include <utility>
#include <vector>

namespace Foam
{

defineTypeNameAndDebug(EikonalSolver1D, 0);
addToRunTimeSelectionTable
(
    conductionSystemSolver,
    EikonalSolver1D,
    dictionary
);


EikonalSolver1D::EikonalSolver1D(const fvMesh&, const dictionary& solverCoeffs)
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
    // Dijkstra wavefront propagation on the conduction graph.
    //
    // Activation time at node j: Tact[j] = Tact[i] + L_ij / c0
    //
    // O((N+E) log N) via a min-heap, versus O(N*E) for Bellman-Ford.
    // Correct for any connected graph with non-negative edge weights.

    const conductionGraph& G = domain.graph();
    scalarField& Tact = domain.activationTime();

    if (gMax(Tact) < 0)
    {
        WarningInFunction
            << "No nodes in the conduction system are activated "
            << "(Tact < 0 for all nodes). Eikonal propagation skipped."
            << endl;
        return;
    }

    const scalar c0 = c0_.value();
    const label  N  = G.nNodes;

    // Build undirected adjacency list: adj[node] = list of (neighbour, edgeIdx).
    // Stored as parallel lists to stay compatible with OpenFOAM's labelList.
    List<DynamicList<label>> adjNode(N);
    List<DynamicList<label>> adjEdge(N);

    forAll(G.edgeNodeA, eI)
    {
        label a = G.edgeNodeA[eI];
        label b = G.edgeNodeB[eI];
        adjNode[a].append(b);  adjEdge[a].append(eI);
        adjNode[b].append(a);  adjEdge[b].append(eI);
    }

    // Min-heap: (Tact_value, nodeID).
    // std::greater gives smallest Tact at the top.
    using Entry = std::pair<scalar, label>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

    // Seed with every pre-activated node.
    forAll(Tact, nodeI)
    {
        if (Tact[nodeI] >= 0)
        {
            pq.push(std::make_pair(Tact[nodeI], nodeI));
        }
    }

    while (!pq.empty())
    {
        scalar t_curr = pq.top().first;
        label  i      = pq.top().second;
        pq.pop();

        // Stale entry: a shorter path to i was already processed.
        if (t_curr > Tact[i]) continue;

        forAll(adjNode[i], nbrI)
        {
            label  j     = adjNode[i][nbrI];
            label  eI    = adjEdge[i][nbrI];
            scalar t_new = Tact[i] + G.edgeLengths[eI] / c0;

            if (Tact[j] < 0 || t_new < Tact[j])
            {
                Tact[j] = t_new;
                pq.push(std::make_pair(t_new, j));
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
