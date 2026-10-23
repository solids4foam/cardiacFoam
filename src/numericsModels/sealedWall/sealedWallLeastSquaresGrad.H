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
    Foam::fv::sealedWallLeastSquaresGrad

Description
    Least-squares gradient of a potential with sealed walls.

    Cells without a sealed or value face use OpenFOAM's face-neighbour least
    squares. A wall cell uses the point-neighbour least squares of
    sealedWallLeastSquaresVectors, in which no sealed face value enters, so
    the gradient takes no data from the walls the current cannot cross.

    The per-patch states come from the field's owner through
    sealedWallPolicy; the field is the one named in the gradient name. With
    no sealed or value patch the result is leastSquares.

    Usage in gradSchemes:
        grad(Vm)                sealedWallLeastSquares;
        grad(Vm.component(0))   sealedWallLeastSquares;

SourceFiles
    sealedWallLeastSquaresGrad.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef Foam_sealedWallLeastSquaresGrad_H
#define Foam_sealedWallLeastSquaresGrad_H

#include "gradScheme.H"

namespace Foam
{
namespace fv
{

class sealedWallLeastSquaresGrad
:
    public fv::gradScheme<scalar>
{
public:

    //- Runtime type information
    TypeName("sealedWallLeastSquares");


    // Constructors

        //- Construct from mesh and Istream
        sealedWallLeastSquaresGrad(const fvMesh& mesh, Istream&)
        :
            gradScheme<scalar>(mesh)
        {}


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //- for optional caching
        virtual tmp<volVectorField> calcGrad
        (
            const volScalarField& vsf,
            const word& name
        ) const;
};


} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
