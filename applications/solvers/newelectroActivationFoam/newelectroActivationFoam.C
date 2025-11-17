/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    newManufacturedFoam

Description
    Solves the reaction-diffusion equation for muscle electrophysiology
    stemming from the mono-domain approach, where the ionic model is run-time
    selectable.

Authors
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/ 
#include "fvCFD.H"
#include "newionicModel.H"
#include "pimpleControl.H"
#include "Field.H"


#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // createFields.H should create: Vm, Iion, conductivity, chi, Cm,
    // electroActivationProperties, ionicModelCoeffs, and the autoPtr<newionicModel> newionicModel
    


    pimpleControl pimple(mesh);

    const label nStates = newionicModel->nEqns();
    Field<Field<scalar>> states
    (
        mesh.nCells(),
        Field<scalar>(nStates, 0.0)
    );

    // Time stepping options
    const Switch solveExplicit
    (
        electroActivationProperties.lookup("solveExplicit")
    );

   
    // External state storage: one state vector per cell.
    // The ionic model decides what each index means.
   

    // For the ionic model to initialise its states based on Vm and position,
    // add a virtual newionicModel::initializeStates(Vm, C, states) later.
    // For now we just start with whatever the model uses internally.

    if (solveExplicit)
    {
        scalar dt = runTime.deltaTValue();
        // Sync dt across processors
        reduce(dt, maxOp<scalar>());
        runTime.setDeltaT(dt);

        while (runTime.loop())
        {
            const scalar currentTime = runTime.value();
            const scalar deltaT      = runTime.deltaTValue();
            const scalar t0          = currentTime - deltaT;

            Info<< "\nTime = " << currentTime << " (explicit)" << nl;

            scalarField& externalStimulusCurrentI = externalStimulusCurrent;
            externalStimulusCurrentI = 0.0;
            if (t0 <= stimulusDuration.value())
            {
                forAll(stimulusCellIDs, cI)
                {
                    const label cellID = stimulusCellIDs[cI];
                    externalStimulusCurrent[cellID] = stimulusIntensity.value();
                }
            }
            externalStimulusCurrent.correctBoundaryConditions();
    
            // 1) Compute ionic current using OLD Vm
            newionicModel->calculateCurrent
            (
                t0,                      // stepStartTime
                deltaT,                  // deltaT
                Vm.internalField(),      // Vm_old
                ionicCurrent,                    // Iion
                states                   // opaque state buffer
            );
            ionicCurrent.correctBoundaryConditions();

            // 2) Solve monodomain PDE for Vm
            solve
            (
                chi*Cm*fvm::ddt(Vm)
              == fvc::laplacian(conductivity, Vm) - chi*Cm*ionicCurrent + externalStimulusCurrent
            );
            Vm.correctBoundaryConditions();

            // 3) Advance ionic states (gating variables) using NEW Vm
            newionicModel->calculateGating
            (
                t0,
                deltaT,
                Vm.internalField(),      // Vm_new
                ionicCurrent,                    // can be used or ignored by the model
                states                   // same opaque state buffer
            );

            // Debug: check Vm and first cell's m,h
            if (mesh.nCells() > 0)
            {
                const label c = 0;
                Info<< "[explicit loop] t=" << currentTime
                    << " Vm(cell0)=" << Vm[c]
                    << "  m=" << states[c][7]
                    << "  h=" << states[c][8]
                    << nl;
            }

            runTime.write();
        }
    }
    else // implicit in Vm, explicit/ODE in ionic model
    {
        while (runTime.loop())
        {
            const scalar currentTime = runTime.value();
            const scalar deltaT      = runTime.deltaTValue();
            const scalar t0          = currentTime - deltaT;

            scalarField& externalStimulusCurrentI = externalStimulusCurrent;
            externalStimulusCurrentI = 0.0;
            if (t0 <= stimulusDuration.value())
            {
                forAll(stimulusCellIDs, cI)
                {
                    const label cellID = stimulusCellIDs[cI];
                    externalStimulusCurrent[cellID] = stimulusIntensity.value();
                }
            }
            externalStimulusCurrent.correctBoundaryConditions();
            Info<< "\nTime = " << currentTime << " (implicit)" << nl;

            // 1) Compute ionic current using OLD Vm
            newionicModel->calculateCurrent
            (
                t0,
                deltaT,
                Vm.internalField(),  // Vm_old
                ionicCurrent,
                states
            );
            ionicCurrent.correctBoundaryConditions();

            // 2) Implicit solve for Vm with PIMPLE
            while (pimple.loop())
            {
                solve
                (
                    chi*Cm*fvm::ddt(Vm)
                 == fvm::laplacian(conductivity, Vm) - chi*Cm*ionicCurrent + externalStimulusCurrent
                );
            }
            Vm.correctBoundaryConditions();

            // 3) Advance ionic model in time (ODE solve) with NEW Vm
            newionicModel->solveODE
            (
                t0,
                deltaT,
                Vm.internalField(),  // Vm_new
                ionicCurrent,
                states
            );
            // Debug: check Vm and first cell's m,h
            if (mesh.nCells() > 0)
            {
                const label c = 0;
                Info<< "[explicit loop] t=" << currentTime
                    << " Vm(cell0)=" << Vm[c]
                    << "  m=" << states[c][7]
                    << "  h=" << states[c][8]
                    << nl;
            }

            
            

#           include "updateActivationTimes.H"

            // --- 5️⃣ Write output fields ---
            if (runTime.writeTime())
            {
                activationVelocity = fvc::grad(
                    1.0 / (activationTime + dimensionedScalar("SMALL", dimTime, SMALL))
                );

                runTime.write();
                runTime.printExecutionTime(Info);
            }
        }


        Info << nl << endl;
        runTime.printExecutionTime(Info);
        Info << "End of simulation" << nl << endl;

        return 0;
    }
}




// ************************************************************************* //


 


