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

#include "ORdGPUCompact.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(ORdGPUCompact, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        ORdGPUCompact,
        dictionary
    );
}

Foam::ORdGPUCompact::ORdGPUCompact
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ORdGPUOpt(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

Foam::ORdGPUCompact::~ORdGPUCompact()
{}

void Foam::ORdGPUCompact::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

// ************************************************************************* //
