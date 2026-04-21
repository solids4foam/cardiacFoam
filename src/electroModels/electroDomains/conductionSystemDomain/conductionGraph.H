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

Class
    Foam::conductionGraph

Description
    Value type holding 1D graph topology for a conduction network:
    nodes, edges, lengths, and conductances. Loaded from a dictionary
    containing a 'conductionEdges' list of (nodeA nodeB length conductance) tuples.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef conductionGraph_H
#define conductionGraph_H

#include "scalarField.H"
#include "labelList.H"
#include "dictionary.H"
#include "DynamicList.H"
#include "boolList.H"
#include <queue>

namespace Foam
{

struct conductionGraph
{
    label       nNodes;
    label       nEdges;
    labelList   edgeNodeA;
    labelList   edgeNodeB;
    scalarField edgeLengths;
    scalarField edgeConductances;

    // Tree topology rooted at node 0 (built during readFromDict).
    // Required by Monodomain1DSolver (Hines algorithm).
    labelList   parentList;    // parent[i] = parent node (-1 for root)
    labelList   orderList;     // BFS order root→leaves (back-substitution pass)
    labelList   reverseOrder;  // BFS order leaves→root (forward-elimination pass)

    conductionGraph()
    :
        nNodes(0),
        nEdges(0)
    {}

    void readFromDict(const dictionary& dict, const bool reportSetup = false)
    {
        const List<scalarList> edgeEntries
        (
            dict.get<List<scalarList>>("conductionEdges")
        );
        nEdges = edgeEntries.size();

        edgeNodeA.setSize(nEdges);
        edgeNodeB.setSize(nEdges);
        edgeLengths.setSize(nEdges);
        edgeConductances.setSize(nEdges);

        label maxNode = 0;
        forAll(edgeEntries, eI)
        {
            const scalarList& e = edgeEntries[eI];
            if (e.size() != 4)
            {
                FatalErrorInFunction
                    << "Each entry in 'conductionEdges' must have 4 values: "
                    << "(nodeA nodeB length conductance). Entry " << eI
                    << " has " << e.size() << " values."
                    << exit(FatalError);
            }
            edgeNodeA[eI]        = static_cast<label>(e[0]);
            edgeNodeB[eI]        = static_cast<label>(e[1]);
            edgeLengths[eI]      = e[2];
            edgeConductances[eI] = e[3];
            maxNode = max(maxNode, max(edgeNodeA[eI], edgeNodeB[eI]));
        }

        nNodes = maxNode + 1;

        if (reportSetup)
        {
            Info<< "Conduction graph: " << nNodes << " nodes, "
                << nEdges << " edges." << endl;
        }

        buildTreeTopology(reportSetup);
    }

private:

    // BFS from node 0 to establish parent–child relationships and the
    // two traversal orderings needed by the Hines two-pass algorithm.
    // Aborts with a fatal error if the graph is disconnected or not a tree.
    void buildTreeTopology(const bool reportSetup)
    {
        if (nNodes == 0) return;

        // Validate: a tree on N nodes has exactly N-1 edges.
        if (nEdges != nNodes - 1)
        {
            FatalErrorInFunction
                << "Conduction graph has " << nNodes << " nodes but "
                << nEdges << " edges. A tree requires exactly N-1 = "
                << (nNodes - 1) << " edges. "
                << "Hines algorithm requires an acyclic (tree) topology."
                << exit(FatalError);
        }

        // Build undirected adjacency list.
        List<DynamicList<label>> adj(nNodes);
        forAll(edgeNodeA, eI)
        {
            adj[edgeNodeA[eI]].append(edgeNodeB[eI]);
            adj[edgeNodeB[eI]].append(edgeNodeA[eI]);
        }

        parentList.setSize(nNodes, -1);
        orderList.setSize(nNodes);
        boolList visited(nNodes, false);

        std::queue<label> q;
        q.push(0);
        visited[0] = true;
        label orderIdx = 0;

        while (!q.empty())
        {
            label node = q.front(); q.pop();
            orderList[orderIdx++] = node;

            forAll(adj[node], nbrI)
            {
                label nbr = adj[node][nbrI];
                if (!visited[nbr])
                {
                    visited[nbr] = true;
                    parentList[nbr] = node;
                    q.push(nbr);
                }
            }
        }

        if (orderIdx != nNodes)
        {
            FatalErrorInFunction
                << "Conduction graph is not connected: BFS from node 0 "
                << "reached " << orderIdx << " of " << nNodes << " nodes."
                << exit(FatalError);
        }

        // reverseOrder: reverse of BFS (leaves→root for forward elimination).
        reverseOrder.setSize(nNodes);
        forAll(orderList, i)
        {
            reverseOrder[i] = orderList[nNodes - 1 - i];
        }

        if (reportSetup)
        {
            Info<< "Conduction graph: tree topology built (root=0)." << endl;
        }
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
