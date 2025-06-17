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
    electroActivationTenTuscherFoam

Description
    This solver is a driver for a single cell electrophysiology model from
    cellml.org. In this case, it is used to drive the Ten Tusscher, Noble,
    Noble, Panfilov, 2004 model.

    This solver extends singleCellElectroActivationFoam to solve a PDE for the
    voltage at the OpenFOAM solver level, while the ionic model is still solved
    within the ionicModel class (assuming a constant voltage per time step).

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interpolationTable.H"
#include "OFstream.H"
#include "ionicModelCellML.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // Create a file for the output
    //OFstream output("output.txt");

    // Loop through time
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Set the external stimulus current for this time
        scalarField& externalStimulusCurrentI = externalStimulusCurrent;
        externalStimulusCurrentI = 0.0;
        if (runTime.value() <= stimulusDuration.value())
        {
            forAll(stimulusCellIDs, cI)
            {
                const label cellID = stimulusCellIDs[cI];

                externalStimulusCurrentI[cellID] += stimulusIntensity.value();
            }
        }

        // Solve the ionic model given the current voltage and calculate the
        // ionic model currents
        scalarField& ionicCurrentI = ionicCurrent;
        ionicCurrentI = 0.0;
        ionicModel.calculateCurrent
        (
            runTime.value() - runTime.deltaTValue(),
            runTime.deltaTValue(),
            Vm.internalField(),
            ionicCurrentI
        );

        // Construct and solve the voltage equation given a known ionic current
        // and external stimulus current
        fvScalarMatrix VmEqn
        (
            beta*Cm*fvm::ddt(Vm)
         ==
            fvm::laplacian(conductivity, Vm)
          - BuenoOrovioScaleFactor*beta*Cm*ionicCurrent // GET RID OF THIS SCALE FACTOR!
          + externalStimulusCurrent
        );

        VmEqn.solve();

        if (runTime.writeTime())
        {
            runTime.write();
            runTime.printExecutionTime(Info);
        }
    }

    Info<< nl << endl;

    runTime.printExecutionTime(Info);

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
