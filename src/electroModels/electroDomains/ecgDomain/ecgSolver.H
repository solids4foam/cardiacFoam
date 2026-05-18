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

\*---------------------------------------------------------------------------*/

#ifndef ecgSolver_H
#define ecgSolver_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

namespace Foam
{

class dictionary;
class ecgDomain;

class ecgSolver
{
public:

    TypeName("ecgSolver");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ecgSolver,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    ecgSolver() = default;

    static autoPtr<ecgSolver> New(const dictionary& dict);

    virtual ~ecgSolver() = default;

    virtual void solve
    (
        ecgDomain& domain,
        scalar t0,
        scalar dt,
        scalarField& values
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
