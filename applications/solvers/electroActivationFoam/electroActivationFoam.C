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
    elactroActivationFoam

Description
    Solves the reaction-diffusion equation for muscle electrophysiology
    stemming from the mono-domain approach.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "BuenoOrovioIonicModel.H"
#include "NashPanfilovActiveTensionModel.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solves the reaction-diffusion equation stemming from the mono-domain "
        "approach for muscle electro-activation"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Calculate the ionic model currents
        scalarField& ionicCurrentI = ionicCurrent;
        ionicCurrentI = 0.0;
        ionicModel.calculateCurrent
        (
            runTime.value() - runTime.deltaTValue(),
            runTime.deltaTValue(),
            Vm.internalField(),
            ionicCurrentI
        );

        // Calculate stimulus current
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

        // Calculate the active tension
        scalarField& activeTensionI = activeTension;
        activeTensionI = 0.0;
        activeTensionModel.calculateTension
        (
            runTime.value() - runTime.deltaTValue(),
            runTime.deltaTValue(),
            Vm.internalField(),
            activeTensionI
        );

        //while (pimple.loop())
        {
            fvScalarMatrix VmEqn
            (
                beta*Cm*fvm::ddt(Vm)
             ==
                fvm::laplacian(conductivity, Vm)
              - BuenoOrovioScaleFactor*beta*Cm*ionicCurrent
              + externalStimulusCurrent
            );

            VmEqn.solve();
        }

#       include "updateActivationTimes.H"

        if (runTime.writeTime())
        {
            runTime.write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
