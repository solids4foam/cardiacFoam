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
    Foam::conductionSystemSolver

Description
    Abstract base class for conduction-system solvers.

    Concrete implementations are
    registered in the runtime selection table under the dictionary key
    "conductionSystemSolver".

SourceFiles
    conductionSystemSolver.C
\*---------------------------------------------------------------------------*/

#ifndef conductionSystemSolver_H
#define conductionSystemSolver_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"

namespace Foam
{

class dictionary;
class ConductionSystemDomain;

class conductionSystemSolver
{
public:

    TypeName("conductionSystemSolver");

    declareRunTimeSelectionTable
    (
        autoPtr,
        conductionSystemSolver,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );

    conductionSystemSolver() = default;

    static autoPtr<conductionSystemSolver> New
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    virtual ~conductionSystemSolver() = default;

    virtual void advance
    (
        ConductionSystemDomain& domain,
        scalar t0,
        scalar dt
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
