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
    exposedPhiETrace

Description
    Heart phiE exposed-face trace from the global phiE.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef exposedPhiETrace_H
#define exposedPhiETrace_H

#include "volFields.H"

namespace Foam
{

scalar interfaceValue
(
    const vector& n,
    const point& Cf,
    const point& CH,
    const point& CB,
    const tensor& aH,
    const tensor& aB,
    const scalar phiH,
    const scalar phiB,
    const vector& gradF
);

label useExposedPhiETrace
(
    volScalarField& heartPhiE,
    const labelUList& faceMap,
    const fvMesh& baseMesh
);

void setExposedPhiETrace
(
    volScalarField& heartPhiE,
    const volScalarField& globalPhiE,
    const labelUList& cellMap,
    const labelUList& faceMap,
    const volTensorField& Ge,
    const volTensorField& sigmaGlobal
);

}

#endif

// ************************************************************************* //
