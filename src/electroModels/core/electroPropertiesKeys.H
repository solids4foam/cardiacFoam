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
    Foam::electroPropertiesKeys

Description
    The keywords each block of electroProperties accepts, checked once when
    the electro model is built so that a misspelt or invented keyword stops
    the run.

    Blocks whose keywords depend on a selected model are not checked here:
    verificationModel, ionicConstantOverrides, constants, initialStates,
    ODESolver, electrodePositions, groundPatches, surfaceCurrentPatches and
    the restitution solver's stimulus.

SourceFiles
    electroPropertiesKeys.C

\*---------------------------------------------------------------------------*/

#ifndef electroPropertiesKeys_H
#define electroPropertiesKeys_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electroPropertiesKeys
{

//- Check electroProperties and its <solverType>Coeffs block.
void check(const dictionary& electroProperties, const word& solverType);

} // End namespace electroPropertiesKeys
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
