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
#include "fieldInit.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    
    ionicModelFDA->initializeFields
    (
        Vm, 
        u1, 
        u2, 
        u3, 
        mesh.C()
    );
    // Initial write of fields
    // Force write of initialized fields at t=0
    Vm.write();
    u1.write();
    u2.write();
    u3.write();
    runTime.write();

Info << "Initial fields written to 0/ directory\n" << endl;

    // Loop through time
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;

        
        // Solve the ionic model for the manufactured solutions with 3 gating variable.
        scalarField& ionicCurrentI = ionicCurrent;
        scalarField& u1I = u1;
        scalarField& u2I = u2;
        scalarField& u3I = u3;


        ionicCurrentI = 0.0;
        u1I = 0.0;
        u2I = 0.0;
        u3I = 0.0;
        ionicModelFDA->calculateCurrent
        (
            runTime.value() - runTime.deltaTValue(),
            runTime.deltaTValue(),
            Vm.internalField(),
            ionicCurrentI,
            u1I,
            u2I,
            u3I

        );
        Vm.correctBoundaryConditions();
        u1.correctBoundaryConditions();
        u2.correctBoundaryConditions();
        u3.correctBoundaryConditions();

        // Construct and solve the voltage equation given a known ionic current
        // and  external stimulus current
        fvScalarMatrix VmEqn
        (
            chi*Cm*fvm::ddt(Vm)
         ==
            fvm::laplacian(conductivity, Vm)
          - chi*ionicCurrent
        );

        VmEqn.solve();
        runTime.write();
        
    }

    Info<< nl << endl;

    runTime.printExecutionTime(Info);

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
