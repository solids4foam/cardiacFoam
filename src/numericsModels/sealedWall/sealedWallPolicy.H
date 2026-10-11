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

Namespace
    Foam::sealedWallPolicy

Description
    Which boundary patches of a potential are sealed walls, and the per-patch
    wall-gradient states the field's owner publishes for the
    sealedWallLeastSquares gradient scheme.

    The states live on the field's mesh as the labelIOList
    wallGradientPatches_<field>, one entry per boundary patch:
    - stock:  the patch takes no part in the wall-cell treatment;
    - sealed: no current crosses the patch; its face values carry no data;
    - value:  the patch value is a boundary datum of the potential.

    The gradient scheme finds the states from its gradient name:
    grad(<field>) and grad(<field>.component(N)) both give <field>.

SourceFiles
    sealedWallPolicy.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef sealedWallPolicy_H
#define sealedWallPolicy_H

#include "volFields.H"

namespace Foam
{
namespace sealedWallPolicy
{

//- Wall-gradient state of a boundary patch
enum patchState : label
{
    stock = 0,
    sealed = 1,
    value = 2
};

//- True for a non-coupled patch on which the potential is sealed
bool isSealed(const fvPatchScalarField& psiP);

//- States of psi: sealed walls sealed, non-coupled fixed-value patches value,
//- every other patch stock; every patch stock when seal is false
labelList wallStates(const volScalarField& psi, const bool seal);

//- True if any state is sealed or value
bool anyWall(const labelUList& states);

//- Field name in a gradient name: grad(<f>) and grad(<f>.component(N)) give <f>
word fieldKey(const word& gradName);

//- Registry name of the states of field key
word statesName(const word& key);

//- Store the states of fieldName on mesh, replacing earlier ones
void publish
(
    const fvMesh& mesh,
    const word& fieldName,
    const labelUList& states
);

//- Store wallStates(psi, seal) for psi on its mesh
void publish(const volScalarField& psi, const bool seal);

//- Published states of field key on mesh; fatal if absent or malformed
const labelList& states(const fvMesh& mesh, const word& key);

} // End namespace sealedWallPolicy
} // End namespace Foam

#endif

// ************************************************************************* //
