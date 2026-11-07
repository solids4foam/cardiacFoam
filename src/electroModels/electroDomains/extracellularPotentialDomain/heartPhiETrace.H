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

Namespace
    Foam::heartPhiETrace

Description
    The heart phiE at the heart-torso wall when there is a torso
    (bathHeartPhiETrace of the <myocardiumSolver>Coeffs block):

    - zeroGradient: the heart phiE has zeroGradient at that wall.
    - global: the heart phiE takes the global phiE on the exposed heart
      faces (exposedPhiETrace); requires sealedHeartBoundary true.

    Read by bidomainSolver (validation), extracellularPotentialDomain (the
    bath sets the trace) and leadFieldConductor (the lead fields use it).

SourceFiles
    heartPhiETrace.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef heartPhiETrace_H
#define heartPhiETrace_H

#include "dictionary.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heartPhiETrace
{

enum class type
{
    zeroGradient,
    global
};

extern const Enum<type> names;

//- bathHeartPhiETrace of the myocardium Coeffs block (required); the single
//  reader of the key. global requires sealedHeartBoundary true.
type readWithTorso(const dictionary& myocardiumCoeffs);

//- The trace of a run: readWithTorso when electroProperties has a torso.
//  Without a torso the key is rejected and the heart wall is the outer,
//  zeroGradient boundary.
type read(const dictionary& myocardiumCoeffs);

} // End namespace heartPhiETrace
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
