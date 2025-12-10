#include "notsure.H"
#include "tmanufacturedFDA.H"   

notsureLoopHandler::notsureLoopHandler
(
    const fvMesh& mesh,
    ionicModel& model
)
:
    mesh_(mesh),
    ionicModel_(model)
{}


// ============================================================================
// Perform one implicit Vm step + ODE evolution using PIMPLE
// ============================================================================
void notsureLoopHandler::notsureLoop
(
    const scalar t0,
    const scalar dt,
    volScalarField& Vm,
    volScalarField& Iion,
    Field<Field<scalar>>& states,
    volScalarField& externalStimulusCurrent,
    const labelList& stimulusCellIDs,
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

    if (t0 <= stimulusDuration)
    {
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
           - chi*Cm*Iion
           + externalStimulusCurrent
        );

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
}
