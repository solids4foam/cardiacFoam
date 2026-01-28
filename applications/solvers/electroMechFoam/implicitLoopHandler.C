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

#include "implicitLoopHandler.H"
#include "tmanufacturedFDA.H"
#include "fvm.H"

implicitLoopHandler::implicitLoopHandler
(
    const fvMesh& mesh,
    ionicModel& model
)
:
    // mesh_(mesh),
    ionicModel_(model)
{}


// ============================================================================
// Perform one implicit Vm step + ODE evolution using PIMPLE
// ============================================================================
void implicitLoopHandler::implicitLoop
(
    const scalar t0,
    const scalar dt,
    volScalarField& Vm,
    volScalarField& Iion,
    Field<Field<scalar>>& states,
    volScalarField& externalStimulusCurrent,
    const List<labelList>& stimulusCellIDsList,
    const List<scalar>& stimulusStartTimes,
    const scalar stimulusIntensity,
    const scalar stimulusDuration,
    const dimensionedScalar& chi,
    const dimensionedScalar& Cm,
    const volTensorField& conductivity,
    pimpleControl& pimple
)
{
    // // 0) Manufactured-specific: store "old" internal states
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).updateStatesOld();
    }

    // 1) External stimulus (internal field only)
    scalarField& externalStimulusCurrentI = externalStimulusCurrent;
    externalStimulusCurrentI = 0.0;

    forAll(stimulusCellIDsList, bI)
    {
        const scalar tStart = stimulusStartTimes[bI];
        if (t0 < tStart || t0 > (tStart + stimulusDuration))
        {
            continue;
        }

        const labelList& stimulusCellIDs = stimulusCellIDsList[bI];
        forAll(stimulusCellIDs, cI)
        {
            const label id = stimulusCellIDs[cI];
            externalStimulusCurrentI[id] = stimulusIntensity;
        }
    }
    externalStimulusCurrent.correctBoundaryConditions();

    // 2) Ionic current using OLD Vm
    ionicModel_.calculateCurrent
    (
        t0,
        dt,
        Vm.internalField(),   // Vm_old
        Iion,
        states
    );
    Iion.correctBoundaryConditions();

    if (ionicModel_.hasManufacturedSolution())
    {
        Iion /= Cm.value();
    }

    // 3) Implicit Vm solve (PIMPLE) - no ODE update inside
    while (pimple.loop())
    {
        solve
        (
            chi*Cm * fvm::ddt(Vm)
          == fvm::laplacian(conductivity, Vm)
           - chi*Cm * Iion
           + externalStimulusCurrent
        );
    }
        Vm.correctBoundaryConditions();

    //4) Manufactured-specific: reset internal states back to OLD
        if (ionicModel_.hasManufacturedSolution())
        {
            refCast<tmanufacturedFDA>(ionicModel_).resetStatesToStatesOld();
        }

    // 5) Advance ionic model in time (ODE solve) with NEW Vm, once per timestep
        ionicModel_.solveODE
        (
            t0,
            dt,
            Vm.internalField(),    // Vm_new
            Iion,
            states
        );

}
