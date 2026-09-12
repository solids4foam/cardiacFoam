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

InNamespace
    Foam

Description
    Shared non-orthogonal corrector loop for the anisotropic conductivity-
    tensor Laplacian solvers. Encodes the single rule "N correctors => N+1
    assemble-and-solve passes", so the deferred non-orthogonal (corrected
    snGrad) contribution is refreshed each pass. One vocabulary for the
    monodomain, bidomain, and bath extracellular-phiE solves.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef nonOrthogonalCorrectorLoop_H
#define nonOrthogonalCorrectorLoop_H

#include "pimpleControl.H"
#include "label.H"

namespace Foam
{

//- Control-driven form: re-assemble and solve while the solutionControl's
//  non-orthogonal corrector is active. Preserves OpenFOAM's own
//  finalNonOrthogonalIter() state. Use where a live pimpleControl exists.
template<class AssembleAndSolve>
inline void correctNonOrthogonalLoop
(
    pimpleControl& pimple,
    AssembleAndSolve&& assembleAndSolve
)
{
    while (pimple.correctNonOrthogonal())
    {
        assembleAndSolve();
    }
}


//- Count-driven form: re-assemble and solve max(nCorr,0)+1 times. Use where
//  no solutionControl loop wraps the solve (e.g. the once-per-step bath phiE).
template<class AssembleAndSolve>
inline void correctNonOrthogonalLoop
(
    const label nCorr,
    AssembleAndSolve&& assembleAndSolve
)
{
    for (label corr = 0; corr <= max(nCorr, label(0)); ++corr)
    {
        assembleAndSolve();
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
