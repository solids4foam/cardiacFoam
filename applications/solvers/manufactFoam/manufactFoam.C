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
    electroActivationFoam

Description
    Solves the reaction-diffusion equation for muscle electrophysiology
    stemming from the mono-domain approach, where the ionic model is run-time
    selectable.

Authors
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ionicModelFDA.H"

#include "pimpleControl.H"
#include "manufacturedFields.H"
#include "manufacturedFDA.H"
#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    pimpleControl pimple(mesh);


    // Initialise fields
    ionicModelFDA->initializeFields
    (
        Vm,
        u1,
        u2,
        u3,
        mesh.C()
    );

    const Switch solveExplicit
    (
        electroActivationProperties.lookup("solveExplicit")
    );

    const int dim = ionicModelFDA->tissue();
    Info << dim << endl;
    const double totalCells = returnReduce(mesh.nCells(), sumOp<int>());
    Info << "Total number of cells: " << totalCells << endl;
    const int N = std::round(Foam::pow(totalCells, 1.0 / dim));

    const scalar cfl = 0.1;
    const scalarField x(mesh.C().component(vector::X));
    const double dx = 1.0 / double(N);
    scalar dt = runTime.deltaTValue();
    int nsteps = int(std::ceil(runTime.endTime().value() / dt));

    if (solveExplicit)
    {
	// Explicit diffusion stability
	dt = cfl * dx*dx / max(conductivity.component(tensor::XX)).value();
	nsteps = int(std::ceil(runTime.endTime().value() / dt));
	dt = runTime.endTime().value() / double(nsteps); // snap to hit Tfinal

	// Sync the dt across processors
	reduce(dt, maxOp<scalar>());

	runTime.setDeltaT(dt);

	while (runTime.loop())
	{
	    refCast<manufacturedFDA>(*ionicModelFDA).updateStatesOld();

	    // --- 1️⃣ Compute Iion for PDE solve using old Vm ---
	    ionicModelFDA->calculateCurrent
	    (
		runTime.value() - runTime.deltaTValue(), // stepStartTime
		runTime.deltaTValue(),                  // deltaT
		Vm.internalField(),                     // input Vm_old
		Iion,                                   // output Iion
		u1, u2, u3                              // output gating vars (not advanced yet)
	    );
	    Iion.correctBoundaryConditions();

            //Pout<< __FILE__<< __LINE__ << endl;

	    // --- 2️⃣ Solve diffusion PDE ---
	    solve
	    (
		chi*Cm*fvm::ddt(Vm) == fvc::laplacian(conductivity, Vm) - chi*Iion
	    );



	    ionicModelFDA->calculateGating(
		runTime.value() - runTime.deltaTValue(), // stepStartTime
		runTime.deltaTValue(),                  // deltaT
		Vm.internalField(),                     // input Vm_new
		Iion,                                   // output Iion (optional update)
		u1, u2, u3                              // update gating variables
	    );

	    //Pout<< __FILE__<< __LINE__ << endl;

	    u1.correctBoundaryConditions();
	    u2.correctBoundaryConditions();
	    u3.correctBoundaryConditions();

	    runTime.write();
	}
	}

    else // solve implicit
	{
	// Info<< "deltaT (before) = " << runTime.deltaTValue() << endl;
	// runTime.setDeltaT(dt);
	// Info<< "deltaT (after) = " << runTime.deltaTValue() << endl;

	while (runTime.loop())
	{
	    Info<< nl << "Time = " << runTime.value() << endl;

	    // Update ionic current explicitly
	    refCast<manufacturedFDA>(*ionicModelFDA).updateStatesOld();
	    Info << "calculate Iion" << endl;
	    // --- 1️⃣ Compute Iion for PDE solve using old Vm ---
	    ionicModelFDA->calculateCurrent
	    (
		runTime.value() - runTime.deltaTValue(), // stepStartTime
		runTime.deltaTValue(),                  // deltaT
		Vm.internalField(),                     // input Vm_old
		Iion,                                   // output Iion
		u1, u2, u3                              // output gating vars (not advanced yet)
	    );
	    Iion.correctBoundaryConditions();

	    // Outer iteration implicit loop
	    while (pimple.loop())
	    {
		// Update Vm
		solve
		(
		    chi*Cm*fvm::ddt(Vm)
		 == fvm::laplacian(conductivity, Vm)
		  - chi*Iion
		);
	    }

	    // Update ionic model explicitly

	    // Update the old-time STATES within the ionic model
	    //// Before solving the ionic model, reset its state to the old time
	    refCast<manufacturedFDA>(*ionicModelFDA).resetStatesToStatesOld();
	    ionicModelFDA->solveODE
	    (
		runTime.value() - runTime.deltaTValue(),
		runTime.deltaTValue(),
		Vm.internalField(),
		Iion,
		u1,
		u2,
		u3
	    );
	    Vm.correctBoundaryConditions();
	    u1.correctBoundaryConditions();
	    u2.correctBoundaryConditions();
	    u3.correctBoundaryConditions();


	}
	}

	scalar Tfinal = runTime.value();


	computeAndPrintErrors
	(
	    Vm.internalField(),
	    u1.internalField(),
	    u2.internalField(),
	    mesh.C().component(vector::X),
	    mesh.C().component(vector::Y),
	    mesh.C().component(vector::Z),
	    Tfinal,
	    dim,            // dimension
	    N,
	    dx,
	    dt,
	    nsteps,
	    solveExplicit //,
	    //""
	);

    Info<< "End" << nl << endl;

    return 0;
}










// ************************************************************************* //

