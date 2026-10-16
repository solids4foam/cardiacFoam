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

Function
    conormalStabilisation

Description
    Term added to the Laplacian of a potential with a tensor face
    conductivity Gf so that its non-orthogonal and cross-diffusion correctors
    converge when the conductivity is oblique to the faces:

        fvm::laplacian(dk, psi) - fvc::laplacian(dk, psi),
        dk = |Gf & n| - n & Gf & n

    The implicit face coefficient of the stabilised operator is the conormal
    magnitude |Gf & n|. The term vanishes at convergence of the correctors,
    and dk is zero on faces whose normal is a principal direction of Gf.

SourceFiles
    conormalStabilisation.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef conormalStabilisation_H
#define conormalStabilisation_H

#include "fvMatrices.H"
#include "surfaceFields.H"

namespace Foam
{

//- Stabilisation term of the Laplacian of psi with face conductivity Gf
tmp<fvScalarMatrix> conormalStabilisation
(
    const surfaceTensorField& Gf,
    const volScalarField& psi
);

}

#endif

// ************************************************************************* //
