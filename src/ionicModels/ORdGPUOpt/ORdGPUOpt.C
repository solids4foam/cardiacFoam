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

#include "ORdGPUOpt.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(ORdGPUOpt, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        ORdGPUOpt,
        dictionary
    );
}

Foam::ORdGPUOpt::ORdGPUOpt
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ORdGPU
    (
        dict,
        num,
        initialDeltaT,
        solveVmWithinODESolver,
        FullEvaluationBackend::optimized
    )
{}

Foam::ORdGPUOpt::~ORdGPUOpt()
{}

void Foam::ORdGPUOpt::solveODE
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
