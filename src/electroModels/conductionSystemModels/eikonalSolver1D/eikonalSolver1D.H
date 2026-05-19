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
    eikonalSolver1D

Description
    Eikonal-equation activation-time solver for the 1D conduction system graph.

SourceFiles
    eikonalSolver1D.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef eikonalSolver1D_H
#define eikonalSolver1D_H

#include "dimensionedScalar.H"
#include "conductionSystemSolver.H"
#include "fvMesh.H"
#include "overrideTypeName.H"
#include "scalarField.H"

namespace Foam
{

class eikonalSolver1D
:
    public conductionSystemSolver
{
    dimensionedScalar c0_;        //- conduction velocity [m/s]

public:

    OverrideTypeName("eikonalSolver1D");

    eikonalSolver1D(const fvMesh&, const dictionary& solverCoeffs);

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
