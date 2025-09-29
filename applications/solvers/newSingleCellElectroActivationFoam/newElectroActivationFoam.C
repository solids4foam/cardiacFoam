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
    newElectroActivationFoam

Description
    TO BE UPDATED SIMAO!!!!! Our new single cell solver!
    

    Solves the reaction-diffusion equation for muscle electrophysiology
    stemming from the mono-domain approach, where the ionic model is run-time
    selectable.

Authors
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ionicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"

    // Create commonly used dimensionSets for convenience
    // Read electoActivationProperties
    Info<< "Reading electroActivationProperties\n" << endl;
    IOdictionary electroActivationProperties
    (
        IOobject
        (
            "electroActivationProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Read ionic model dict
    const dictionary& ionicModelCoeffs =
        electroActivationProperties.subDict("ionicModelCoeffs");

    // Create ionicModelCellML object
    autoPtr<ionicModel> ionicModel =
        ionicModel::New
        (
            ionicModelCoeffs,
            1, // number of integration points
            runTime.deltaTValue(),
            true // solve Vm equation within ODE system
        );

    // Loop through time
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;


        // Solve the ionic model given the current voltage and calculate the
        // ionic model currents

        scalarField dummyVmField(1, 0);
        scalarField dummyIonicCurrentField(1, 0);
        ionicModel->calculateCurrent
        (
            runTime.value() - runTime.deltaTValue(),
            runTime.deltaTValue(),
            dummyVmField,
            dummyIonicCurrentField
        );
    }

    Info<< nl << endl;

    runTime.printExecutionTime(Info);

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
