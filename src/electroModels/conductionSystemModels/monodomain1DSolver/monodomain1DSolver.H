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
    Implicit monodomain cable-equation solver for the conduction system
    graph.  Uses the same operator-splitting scheme as the 3-D
    MyocardiumDomain:

        1.  Solve the ionic-model ODE  (reaction step)
        2.  Solve the implicit cable equation  (diffusion step)

    The cable equation on the graph is discretised in exactly the same way
    that fvm::laplacian + fvm::ddt would produce on a 1-D finite-volume
    mesh.  Because the graph is small (typically < 1000 nodes), the
    resulting dense linear system is solved with LU decomposition.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    Monodomain1DSolver.C
\*---------------------------------------------------------------------------*/

#ifndef MonoDomain1DSolver_H
#define MonoDomain1DSolver_H

#include "graphConductionSystemSolver.H"
#include "overrideTypeName.H"
#include "scalarField.H"
#include "scalarMatrices.H"

namespace Foam
{

class Monodomain1DSolver
:
    public GraphConductionSystemSolver
{
    mutable scalarField appliedCurrentBuffer_;

public:

    OverrideTypeName("monodomain1DSolver");

    explicit Monodomain1DSolver(const dictionary&)
    {}

    virtual void advance
    (
        GraphConductionSystemDomain& domain,
        scalar t0,
        scalar dt
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
