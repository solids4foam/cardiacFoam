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

#include "probedOperator.H"
#include "syncTools.H"
#include "globalMeshData.H"
#include "PstreamBuffers.H"
#include "DynamicList.H"
#include "Map.H"

namespace Foam
{

namespace
{

// Rows of all processors, concatenated in processor order (global order)
labelListList gatherRows(const labelListList& local)
{
    List<labelListList> all(Pstream::nProcs());
    all[Pstream::myProcNo()] = local;
    Pstream::gatherList(all);

    labelListList rows;
    if (Pstream::master())
    {
        label n = 0;
        for (const labelListList& r : all) n += r.size();
        rows.setSize(n);
        n = 0;
        for (labelListList& r : all)
        {
            for (labelList& row : r) rows[n++].transfer(row);
        }
    }
    return rows;
}

} // End anonymous namespace


labelListList probedOperator::faceGraph() const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList cellGlobal
    (
        identity(mesh_.nCells(), globalCells_.localStart())
    );
    labelList nbrGlobal;
    syncTools::swapBoundaryCellList(mesh_, cellGlobal, nbrGlobal);

    List<DynamicList<label>> rows(mesh_.nCells());
    forAll(nei, faceI)
    {
        rows[own[faceI]].append(cellGlobal[nei[faceI]]);
        rows[nei[faceI]].append(cellGlobal[own[faceI]]);
    }
    for (const polyPatch& pp : patches)
    {
        if (!pp.coupled()) continue;

        const labelUList& faceCells = pp.faceCells();
        forAll(pp, i)
        {
            const label bFaceI = pp.start() + i - mesh_.nInternalFaces();
            rows[faceCells[i]].append(nbrGlobal[bFaceI]);
        }
    }

    labelListList graph(rows.size());
    forAll(rows, cellI)
    {
        graph[cellI].transfer(rows[cellI]);
    }
    return graph;
}


labelListList probedOperator::pointGraph() const
{
    // Global point labels, with points shared by processors merged
    labelList pointToGlobal;
    labelList uniquePoints;
    const autoPtr<globalIndex> globalPointsPtr
    (
        mesh_.globalData().mergePoints(pointToGlobal, uniquePoints)
    );

    // Each processor's cells as rows of global point labels
    const labelListList& cellPoints = mesh_.cellPoints();
    labelListList cellGlobalPoints(mesh_.nCells());
    forAll(cellPoints, cellI)
    {
        cellGlobalPoints[cellI] =
            labelList(UIndirectList<label>(pointToGlobal, cellPoints[cellI]));
    }
    const labelListList allCellPoints(gatherRows(cellGlobalPoints));

    // On the master: cells of each point, then cells sharing a point
    labelListList graph;
    if (Pstream::master())
    {
        const label nPoints = globalPointsPtr->totalSize();
        List<DynamicList<label>> pointCells(nPoints);
        forAll(allCellPoints, cellI)
        {
            for (const label pointI : allCellPoints[cellI])
            {
                pointCells[pointI].append(cellI);
            }
        }

        graph.setSize(allCellPoints.size());
        labelList mark(allCellPoints.size(), -1);
        forAll(allCellPoints, cellI)
        {
            DynamicList<label> row;
            mark[cellI] = cellI;
            for (const label pointI : allCellPoints[cellI])
            {
                for (const label otherI : pointCells[pointI])
                {
                    if (mark[otherI] != cellI)
                    {
                        mark[otherI] = cellI;
                        row.append(otherI);
                    }
                }
            }
            graph[cellI].transfer(row);
        }
    }
    return graph;
}


labelListList probedOperator::composePattern
(
    const List<labelListList>& graphs,
    const label nCells
)
{
    labelListList pattern(nCells);
    labelList mark(nCells, -1);
    DynamicList<label> row;
    DynamicList<label> next;

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        row.clear();
        row.append(cellI);
        mark[cellI] = cellI;

        // Every hop keeps the cells reached so far (each hop includes self)
        for (const labelListList& graph : graphs)
        {
            next.clear();
            for (const label k : row)
            {
                for (const label j : graph[k])
                {
                    if (mark[j] != cellI)
                    {
                        mark[j] = cellI;
                        next.append(j);
                    }
                }
            }
            row.append(next);
        }

        pattern[cellI] = row;
    }
    return pattern;
}


labelList probedOperator::colourColumns
(
    const labelListList& pattern,
    const label nCells,
    label& nColours
)
{
    // Rows of each column
    labelList nRows(nCells, 0);
    for (const labelList& row : pattern)
    {
        for (const label j : row) ++nRows[j];
    }
    labelListList columnRows(nCells);
    forAll(columnRows, j)
    {
        columnRows[j].setSize(nRows[j]);
        nRows[j] = 0;
    }
    forAll(pattern, i)
    {
        for (const label j : pattern[i])
        {
            columnRows[j][nRows[j]++] = i;
        }
    }

    // Smallest colour not used by a column sharing a row with j
    labelList colour(nCells, -1);
    DynamicList<label> forbidden;
    nColours = 0;
    forAll(colour, j)
    {
        for (const label i : columnRows[j])
        {
            for (const label k : pattern[i])
            {
                if (colour[k] >= 0)
                {
                    while (forbidden.size() <= colour[k])
                    {
                        forbidden.append(-1);
                    }
                    forbidden[colour[k]] = j;
                }
            }
        }

        label c = 0;
        while (c < forbidden.size() && forbidden[c] == j) ++c;
        if (c == forbidden.size())
        {
            forbidden.append(-1);
        }
        colour[j] = c;
        nColours = max(nColours, c + 1);
    }
    return colour;
}


void probedOperator::buildExchange()
{
    const label nProcs = Pstream::nProcs();

    slots_.setSize(columns_.size());
    List<Map<label>> slotOf(nProcs);
    List<DynamicList<label>> requests(nProcs);
    List<DynamicList<label>> recvSlots(nProcs);
    nRemote_ = 0;

    forAll(columns_, i)
    {
        slots_[i].setSize(columns_[i].size(), -1);
        forAll(columns_[i], k)
        {
            const label col = columns_[i][k];
            if (globalCells_.isLocal(col)) continue;

            const label procI = globalCells_.whichProcID(col);
            auto iter = slotOf[procI].cfind(col);
            if (iter.good())
            {
                slots_[i][k] = iter.val();
            }
            else
            {
                slotOf[procI].insert(col, nRemote_);
                requests[procI].append(globalCells_.toLocal(procI, col));
                recvSlots[procI].append(nRemote_);
                slots_[i][k] = nRemote_++;
            }
        }
    }

    recvSlots_.setSize(nProcs);
    forAll(recvSlots, procI)
    {
        recvSlots_[procI].transfer(recvSlots[procI]);
    }

    sendCells_.setSize(nProcs);
    if (Pstream::parRun())
    {
        PstreamBuffers pBufs;
        forAll(requests, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                UOPstream os(procI, pBufs);
                os << labelList(requests[procI]);
            }
        }
        pBufs.finishedSends();
        forAll(sendCells_, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                UIPstream is(procI, pBufs);
                is >> sendCells_[procI];
            }
        }
    }
}


scalarField probedOperator::remoteValues(const scalarField& x) const
{
    scalarField remote(nRemote_, Zero);
    if (!Pstream::parRun())
    {
        return remote;
    }

    PstreamBuffers pBufs;
    forAll(sendCells_, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            UOPstream os(procI, pBufs);
            os << scalarField(UIndirectList<scalar>(x, sendCells_[procI]));
        }
    }
    pBufs.finishedSends();
    forAll(recvSlots_, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            UIPstream is(procI, pBufs);
            const scalarField values(is);
            forAll(values, k)
            {
                remote[recvSlots_[procI][k]] = values[k];
            }
        }
    }
    return remote;
}


probedOperator::probedOperator
(
    const fvMesh& mesh,
    const List<hop>& hops,
    const operatorType& op
)
:
    mesh_(mesh),
    globalCells_(mesh.nCells()),
    nColours_(0),
    columns_(),
    values_(),
    slots_(),
    recvSlots_(),
    sendCells_(),
    nRemote_(0)
{
    const label nCells = globalCells_.totalSize();

    // Hop graphs over all cells, on the master
    const labelListList allFace(gatherRows(faceGraph()));
    List<labelListList> graphs(hops.size());
    bool hasPoint = false;
    forAll(hops, hopI)
    {
        hasPoint = hasPoint || hops[hopI] == hop::point;
    }
    labelListList allPoint;
    if (hasPoint)
    {
        allPoint = pointGraph();
    }
    if (Pstream::master())
    {
        forAll(hops, hopI)
        {
            graphs[hopI] = hops[hopI] == hop::face ? allFace : allPoint;
        }
    }

    // Pattern and colours on the master; each processor receives its rows
    // with the colour of every column, and the colours of its cells
    labelListList rowColours;
    labelList cellColour;
    {
        labelListList pattern;
        labelList colour;
        if (Pstream::master())
        {
            pattern = composePattern(graphs, nCells);
            colour = colourColumns(pattern, nCells, nColours_);
        }
        Pstream::broadcast(nColours_);

        PstreamBuffers pBufs;
        if (Pstream::master())
        {
            for (label procI = 0; procI < Pstream::nProcs(); ++procI)
            {
                const label start = globalCells_.localStart(procI);
                const label size = globalCells_.localSize(procI);
                labelListList procRows(size);
                labelListList procColours(size);
                forAll(procRows, i)
                {
                    procRows[i] = pattern[start + i];
                    procColours[i] =
                        labelList(UIndirectList<label>(colour, procRows[i]));
                }
                const labelList procCellColour
                (
                    SubList<label>(colour, size, start)
                );

                if (procI == Pstream::myProcNo())
                {
                    columns_.transfer(procRows);
                    rowColours.transfer(procColours);
                    cellColour = procCellColour;
                }
                else
                {
                    UOPstream os(procI, pBufs);
                    os << procRows << procColours << procCellColour;
                }
            }
        }
        pBufs.finishedScatters();
        if (!Pstream::master())
        {
            UIPstream is(Pstream::masterNo(), pBufs);
            is >> columns_ >> rowColours >> cellColour;
        }
    }

    // Entries of each colour: (row, position in row)
    List<DynamicList<labelPair>> colourEntries(nColours_);
    values_.setSize(columns_.size());
    forAll(rowColours, i)
    {
        values_[i].setSize(rowColours[i].size(), Zero);
        forAll(rowColours[i], k)
        {
            colourEntries[rowColours[i][k]].append(labelPair(i, k));
        }
    }

    // One probe per colour
    scalarField x(mesh_.nCells());
    for (label c = 0; c < nColours_; ++c)
    {
        forAll(x, cellI)
        {
            x[cellI] = cellColour[cellI] == c ? 1.0 : 0.0;
        }
        const tmp<scalarField> tCx(op(x));
        const scalarField& Cx = tCx();
        for (const labelPair& entry : colourEntries[c])
        {
            values_[entry.first()][entry.second()] = Cx[entry.first()];
        }
    }

    buildExchange();
}


label probedOperator::nNonZeros() const
{
    label n = 0;
    for (const labelList& row : columns_) n += row.size();
    return returnReduce(n, sumOp<label>());
}


tmp<scalarField> probedOperator::apply(const scalarField& x) const
{
    const scalarField remote(remoteValues(x));
    tmp<scalarField> tCx(new scalarField(columns_.size(), Zero));
    scalarField& Cx = tCx.ref();

    forAll(columns_, i)
    {
        forAll(columns_[i], k)
        {
            const label slot = slots_[i][k];
            const scalar xj =
                slot >= 0
              ? remote[slot]
              : x[globalCells_.toLocal(columns_[i][k])];
            Cx[i] += values_[i][k]*xj;
        }
    }
    return tCx;
}


tmp<scalarField> probedOperator::applyAbs(const scalarField& x) const
{
    const scalarField remote(remoteValues(x));
    tmp<scalarField> tCx(new scalarField(columns_.size(), Zero));
    scalarField& Cx = tCx.ref();

    forAll(columns_, i)
    {
        forAll(columns_[i], k)
        {
            const label slot = slots_[i][k];
            const scalar xj =
                slot >= 0
              ? remote[slot]
              : x[globalCells_.toLocal(columns_[i][k])];
            Cx[i] += mag(values_[i][k])*mag(xj);
        }
    }
    return tCx;
}


namespace
{

// Add the off-processor contributions to the processors owning them
void sendContributions
(
    const scalarField& remote,
    const labelListList& recvSlots,
    const labelListList& sendCells,
    scalarField& result
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    PstreamBuffers pBufs;
    forAll(recvSlots, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            UOPstream os(procI, pBufs);
            os << scalarField(UIndirectList<scalar>(remote, recvSlots[procI]));
        }
    }
    pBufs.finishedSends();
    forAll(sendCells, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            UIPstream is(procI, pBufs);
            const scalarField values(is);
            forAll(values, k)
            {
                result[sendCells[procI][k]] += values[k];
            }
        }
    }
}

} // End anonymous namespace


tmp<scalarField> probedOperator::applyTranspose(const scalarField& z) const
{
    tmp<scalarField> tCtz(new scalarField(columns_.size(), Zero));
    scalarField& Ctz = tCtz.ref();
    scalarField remote(nRemote_, Zero);

    forAll(columns_, i)
    {
        forAll(columns_[i], k)
        {
            const scalar v = values_[i][k]*z[i];
            const label slot = slots_[i][k];
            if (slot >= 0)
            {
                remote[slot] += v;
            }
            else
            {
                Ctz[globalCells_.toLocal(columns_[i][k])] += v;
            }
        }
    }

    sendContributions(remote, recvSlots_, sendCells_, Ctz);
    return tCtz;
}


tmp<scalarField> probedOperator::applyTransposeAbs(const scalarField& z) const
{
    tmp<scalarField> tCtz(new scalarField(columns_.size(), Zero));
    scalarField& Ctz = tCtz.ref();
    scalarField remote(nRemote_, Zero);

    forAll(columns_, i)
    {
        forAll(columns_[i], k)
        {
            const scalar v = mag(values_[i][k])*mag(z[i]);
            const label slot = slots_[i][k];
            if (slot >= 0)
            {
                remote[slot] += v;
            }
            else
            {
                Ctz[globalCells_.toLocal(columns_[i][k])] += v;
            }
        }
    }

    sendContributions(remote, recvSlots_, sendCells_, Ctz);
    return tCtz;
}

} // End namespace Foam

// ************************************************************************* //
