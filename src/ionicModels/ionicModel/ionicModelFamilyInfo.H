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

#ifndef ionicModelFamilyInfo_H
#define ionicModelFamilyInfo_H

#include "label.H"
#include "scalar.H"
#include "ionicModelIO.H"

namespace Foam
{

struct ionicModelFamilyInfo
{
    label nConstants;
    label nStates;
    label nAlgebraics;
    const char* const* constantNames;
    const char* const* stateNames;
    const char* const* algebraicNames;
    label voltageStateI;
    scalar timeScale;
    scalar vmInputScale;
    scalar vmInputOffset;
    ionicModelIO::VmTransform vmTransform;
};

} // End namespace Foam

#endif

// ************************************************************************* //
