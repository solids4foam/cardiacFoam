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

#ifndef bathECGSolver_H
#define bathECGSolver_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class dictionary;
class bathDomain;

class bathECGSolver
{
public:

    TypeName("bathECGSolver");

    declareRunTimeSelectionTable
    (
        autoPtr,
        bathECGSolver,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    bathECGSolver() = default;

    static autoPtr<bathECGSolver> New(const dictionary& dict);

    virtual ~bathECGSolver() = default;

    virtual void solve
    (
        bathDomain& domain,
        scalar t0,
        scalar dt
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
