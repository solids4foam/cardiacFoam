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
    Foam::sealedWallLeastSquaresVectors

Description
    Least-squares vectors of the wall cells of one field, cached on the mesh
    under the field's name.

    A wall cell has a face on a sealed or value patch. Its stencil is
    OpenFOAM's point-neighbour stencil (centredCPCCellToCellStencilObject):
    every cell sharing a vertex with it and every non-coupled boundary face
    sharing a vertex with it. Faces of sealed patches are dropped; every
    other entry is a row with weight 1/|d|^2, as in pointCellsLeastSquares.

SourceFiles
    sealedWallLeastSquaresVectors.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef sealedWallLeastSquaresVectors_H
#define sealedWallLeastSquaresVectors_H

#include "MeshObject.H"
#include "fvMesh.H"

namespace Foam
{

class sealedWallLeastSquaresVectors
:
    public MeshObject<fvMesh, MoveableMeshObject, sealedWallLeastSquaresVectors>
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            MoveableMeshObject,
            sealedWallLeastSquaresVectors
        > MeshObject_type;


    // Private Data

        //- Per-patch states the vectors are built from
        const labelList states_;

        //- Wall cells
        labelList wallCells_;

        //- Least-squares vectors of each wall cell, one per stencil entry
        List<List<vector>> vectors_;


    // Private Member Functions

        //- Build the wall cells and their vectors
        void calcVectors();


public:

    // Declare name of the class and its debug switch
    TypeName("sealedWallLeastSquaresVectors");


    // Constructors

        //- Construct with registration name, mesh and per-patch states
        sealedWallLeastSquaresVectors
        (
            const word& name,
            const fvMesh& mesh,
            const labelUList& states
        );


    //- Destructor
    virtual ~sealedWallLeastSquaresVectors() = default;


    // Member Functions

        //- Per-patch states the vectors are built from
        const labelList& states() const
        {
            return states_;
        }

        //- Wall cells
        const labelList& wallCells() const
        {
            return wallCells_;
        }

        //- Least-squares vectors of each wall cell
        const List<List<vector>>& vectors() const
        {
            return vectors_;
        }

        //- Rebuild the vectors when the mesh moves
        virtual bool movePoints();
};


} // End namespace Foam

#endif

// ************************************************************************* //
