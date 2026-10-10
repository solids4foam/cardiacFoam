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
    Foam::interfaceGradient

Description
    Cell gradient of a potential whose normal derivative jumps across the
    heart-bath interface (the faces between cells with and without sigmaI).
    Cells next to the interface use least squares on their own side, with the
    face value from two-sided flux continuity (interfaceValue) on interface
    faces; the face values use the tangential part of gradPrev. Every other
    cell keeps fvc::grad.

SourceFiles
    interfaceGradient.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceGradient_H
#define interfaceGradient_H

#include "volFields.H"

namespace Foam
{

tmp<volVectorField> interfaceGradient
(
    const volScalarField& phi,
    const volTensorField& sigmaE,
    const volTensorField& sigmaI,
    const volVectorField& gradPrev
);

}

#endif

// ************************************************************************* //
