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
    Foam::EikonalSolver1D

Description
    Eikonal-equation activation-time solver for the 1D conduction system
    graph (Purkinje network).

    Computes nodal activation times by propagating a wavefront at a
    prescribed conduction velocity c0 through the graph adjacency list.
    No ionic model is required — the single parameter c0 [m/s] is read
    from the solver coefficients sub-dictionary.

    This solver is the 1D counterpart of EikonalSolver (3D myocardium)
    and is registered in the GraphConductionSystemSolver runtime table
    under the key "eikonalSolver".

Usage
    In electroProperties:
    \verbatim
    conductionSystemSolver  eikonalSolver;
    eikonalSolverCoeffs
    {
        c0  [0 1 -1 0 0 0 0] 0.8;  // conduction velocity [m/s]
    }
    \endverbatim

SourceFiles
    eikonalSolver1D.C
\*---------------------------------------------------------------------------*/

#ifndef EikonalSolver1D_H
#define EikonalSolver1D_H

#include "graphConductionSystemSolver.H"
#include "overrideTypeName.H"
#include "scalarField.H"

namespace Foam
{

class EikonalSolver1D
:
    public GraphConductionSystemSolver
{
    scalar c0_;        //- conduction velocity [m/s]

public:

    OverrideTypeName("eikonalSolver");

    explicit EikonalSolver1D(const dictionary& solverCoeffs);

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
