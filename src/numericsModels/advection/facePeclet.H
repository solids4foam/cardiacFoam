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
    facePeclet

Description
    Face Peclet number of an advection-diffusion operator with face
    advective flux phi and tensor face diffusivity Gf:
        Pe = |phi|/(gamma |Sf| nonOrthDeltaCoeffs),  gamma = n & Gf & n
    the ratio of the advective to the stock implicit diffusive face coefficient
    n & Gf & n; conormal stabilisation raises the implicit coefficient above it.
    Above 2 a centred or higher-order convection scheme is no longer bounded
    by the diffusion, and the solution is under-resolved on that face.
    reportFacePeclet logs its maximum, mean and the share of internal faces
    above 2.

SourceFiles
    facePeclet.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef facePeclet_H
#define facePeclet_H

#include "surfaceFields.H"

namespace Foam
{

//- Peclet number of the internal faces; zero on boundary faces
tmp<surfaceScalarField> facePeclet
(
    const surfaceScalarField& phi,
    const surfaceTensorField& Gf
);


//- Log the maximum, the mean and the share of internal faces above 2
void reportFacePeclet(const surfaceScalarField& Pe, const word& fieldName);

}

#endif

// ************************************************************************* //
