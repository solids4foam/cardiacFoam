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

#ifndef ECGSolver_H
#define ECGSolver_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

namespace Foam
{

class dictionary;
class ECGDomain;

class ECGSolver
{
public:

    TypeName("ECGSolver");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ECGSolver,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    ECGSolver() = default;

    static autoPtr<ECGSolver> New(const dictionary& dict);

    virtual ~ECGSolver() = default;

    virtual void compute
    (
        const ECGDomain& domain,
        scalarField& values
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
