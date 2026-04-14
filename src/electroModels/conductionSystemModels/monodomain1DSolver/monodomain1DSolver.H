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
    Foam::Monodomain1DSolver

Description
    Implicit monodomain cable-equation solver for a tree conduction-system
    graph. Uses the same reaction/diffusion splitting as the 3-D myocardium:

        1.  Solve the ionic-model ODE  (reaction step)
        2.  Solve the implicit cable equation  (diffusion step)

    The cable equation is discretised as a 1-D finite-volume graph problem.
    Each graph node owns a control length assembled from half of each incident
    edge length. Edge diffusion uses sigma/L flux coefficients and is solved
    with a tree Hines elimination pass.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    Monodomain1DSolver.C
\*---------------------------------------------------------------------------*/

#ifndef MonoDomain1DSolver_H
#define MonoDomain1DSolver_H

#include "conductionSystemSolver.H"
#include "fvMesh.H"
#include "overrideTypeName.H"
#include "scalarField.H"
#include "scalarMatrices.H"

namespace Foam
{

class Monodomain1DSolver
:
    public conductionSystemSolver
{
    mutable scalarField appliedCurrentBuffer_;

public:

    OverrideTypeName("monodomain1DSolver");

    Monodomain1DSolver(const fvMesh&, const dictionary&)
    {}

    virtual void advance
    (
        ConductionSystemDomain& domain,
        scalar t0,
        scalar dt
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
