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
    monodomain1DSolver

Description
    Implicit monodomain cable-equation solver for a tree conduction-system graph.

Author
    Simao Nieto de Castro. All rights reserved.

SourceFiles
    monodomain1DSolver.C
\*---------------------------------------------------------------------------*/

#ifndef monodomain1DSolver_H
#define monodomain1DSolver_H

#include "conductionSystemSolver.H"
#include "fvMesh.H"
#include "overrideTypeName.H"
#include "scalarField.H"
#include "scalarMatrices.H"

namespace Foam
{

class monodomain1DSolver
:
    public conductionSystemSolver
{
    mutable scalarField appliedCurrentBuffer_;

public:

    OverrideTypeName("monodomain1DSolver");

    monodomain1DSolver(const fvMesh&, const dictionary&)
    {}

    virtual void advance
    (
        conductionSystemDomain& domain,
        scalar t0,
        scalar dt
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
