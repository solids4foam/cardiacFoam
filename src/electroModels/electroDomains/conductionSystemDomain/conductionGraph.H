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
    containing an 'edges' list of (nodeA nodeB length conductance) tuples.

Author
    Simao Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef conductionGraph_H
#define conductionGraph_H

#include "scalarField.H"
#include "labelList.H"
#include "dictionary.H"

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

    conductionGraph()
    :
        nNodes(0),
        nEdges(0)
    {}

    void readFromDict(const dictionary& dict)
    {
        const List<scalarList> edgeEntries
        (
            dict.get<List<scalarList>>("edges")
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
                    << "Each entry in 'edges' must have 4 values: "
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

        Info<< "Conduction graph: " << nNodes << " nodes, "
            << nEdges << " edges." << endl;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
