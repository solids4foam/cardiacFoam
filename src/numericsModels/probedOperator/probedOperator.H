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
    Foam::probedOperator

Description
    A linear cell-to-cell operator C, given as a function x -> C x on the
    cells of a mesh, stored as a sparse matrix by graph-coloured probing.

    The sparsity pattern is a composition of cell-to-cell hops, outermost
    first:
    - face:  the cell and its face neighbours, across internal faces and
             every coupled patch;
    - point: the cell and the cells that share a point with it.
    Columns that share no row get the same colour (distance-2 colouring of
    the pattern), and C x_c with x_c the indicator of colour c gives every
    entry of that colour. An operator whose stencil exceeds the pattern, or
    that is not linear, is not reproduced; reproductionError() measures it.

    Rows are stored on their owning processor with global column ids. The
    pattern and the colouring are built on the master, in global cell order,
    so the master holds nCells times the row length labels.

    apply() and applyTranspose() exchange off-processor column values and
    contributions with the processors that own them.

SourceFiles
    probedOperator.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef probedOperator_H
#define probedOperator_H

#include "fvMesh.H"
#include "globalIndex.H"
#include <functional>

namespace Foam
{

class probedOperator
{
public:

    //- One hop of the sparsity pattern
    enum class hop { face, point };

    //- The operator: x (one value per local cell) -> C x
    typedef std::function<tmp<scalarField>(const scalarField&)> operatorType;


private:

    const fvMesh& mesh_;

    globalIndex globalCells_;

    label nColours_;

    // Row i: global column ids and values
    labelListList columns_;
    List<scalarList> values_;

    // Off-processor columns: compact slot of each entry (-1 for a local
    // column), and per processor the slots it owns and the local cells it
    // requests from this processor
    labelListList slots_;
    labelListList recvSlots_;
    labelListList sendCells_;
    label nRemote_;

    // Cell graphs in global ids, rows in local cell order
    labelListList faceGraph() const;
    labelListList pointGraph() const;

    // Pattern rows of the global cells, on the master
    static labelListList composePattern
    (
        const List<labelListList>& graphs,
        const label nCells
    );

    // Greedy distance-2 colouring of the pattern columns, global order
    static labelList colourColumns
    (
        const labelListList& pattern,
        const label nCells,
        label& nColours
    );

    void buildExchange();

    // Values of x at the off-processor columns, by compact slot
    scalarField remoteValues(const scalarField& x) const;


public:

    //- Probe op on mesh with the pattern given by hops, outermost first
    probedOperator
    (
        const fvMesh& mesh,
        const List<hop>& hops,
        const operatorType& op
    );

    probedOperator(const probedOperator&) = delete;
    void operator=(const probedOperator&) = delete;


    //- Number of colours (= probes)
    label nColours() const
    {
        return nColours_;
    }

    //- Global number of stored entries
    label nNonZeros() const;

    //- C x
    tmp<scalarField> apply(const scalarField& x) const;

    //- C^T z
    tmp<scalarField> applyTranspose(const scalarField& z) const;

    //- |C| |x|, entry by entry
    tmp<scalarField> applyAbs(const scalarField& x) const;

    //- |C|^T |z|, entry by entry
    tmp<scalarField> applyTransposeAbs(const scalarField& z) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
