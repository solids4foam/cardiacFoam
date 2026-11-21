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
    Foam::torsoCoupling

Description
    How the potential of the region outside the heart is obtained, from the
    top level of electroProperties:

        torso { cellZones (...); conductivityField <name>; ... }
        torsoCoupling bath | leadField;

    - bath: the extracellularPotentialDomain solves phiE on the heart and the
      torso every step, coupled both ways; bidomainSolver only.
    - leadField: torsoECG uses the lead fields that buildLeadFields
      precomputes for the torso; one-way; monodomainSolver or
      bidomainSolver.

    torso and torsoCoupling are given together or not at all.

SourceFiles
    torsoCoupling.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef torsoCoupling_H
#define torsoCoupling_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace torsoCoupling
{

enum class type
{
    none,
    bath,
    leadField
};

//- The torsoCoupling of the top-level electroProperties: none when it has
//  neither torso nor torsoCoupling. Fatal when only one of them is given or
//  the value is not bath or leadField.
type read(const dictionary& electroProperties);

//- read() of the electroProperties that holds myocardiumCoeffs, the
//  <myocardiumSolver>Coeffs block; fatal when myocardiumCoeffs has no
//  enclosing dictionary
type readFromCoeffs(const dictionary& myocardiumCoeffs);

//- The keyword of a coupling (none, bath, leadField)
word name(const type coupling);

} // End namespace torsoCoupling
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
