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
    conductivityFieldIO

Description
    Strict conductivity source resolution shared by the myocardium solvers.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef conductivityFieldIO_H
#define conductivityFieldIO_H

#include "dictionary.H"
#include "fvMeshSubset.H"
#include "volFields.H"

namespace Foam
{

struct conductivityFieldSpec
{
    word fieldName;
    word dictionaryEntry;
    word internalName;
};


tmp<volTensorField> readConductivityField
(
    const fvMesh& solverMesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const dictionary& coefficients,
    const conductivityFieldSpec& spec
);

} // End namespace Foam

#endif

// ************************************************************************* //
