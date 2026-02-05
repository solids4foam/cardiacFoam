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
    singleCellElectroActivationFoam

Description
    Single-cell electrophysiology solver based on OpenFOAM.

    This solver advances a single cardiac cell model in time using a
    run-time selectable ionic model. In contrast to electroActivationFoam,
    no spatial voltage PDE is solved; instead, the transmembrane voltage
    Vm is integrated directly as part of the ionic model ODE system.

    The solver is primarily intended for:
      - single-cell action potential simulations,
      - testing and validation of ionic models,
      - debugging and development of electrophysiology model components.

    External stimulation is prescribed via a stimulus protocol dictionary.
    Output is written as time series of states, algebraic variables and rates
    using ionicModelIO.

    This solver represents a special case of electroActivationFoam. Future
    versions of electroActivationFoam are expected to support single-cell
    operation directly, at which point this solver may be deprecated.

Author
    Simao Nieto de Castro, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "activeTensionModel.H"
#include "activeTensionIO.H"

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createFields.H"

    // Create ionic model
    autoPtr<ionicModel> ionicModel =
        ionicModel::New
        (
            solutionVariablesMemory,
            1, // one integration point
            runTime.deltaTValue(),
            true // solve ODE system only
        );

    autoPtr<activeTensionModel> activeTensionModel =
        activeTensionModel::New
        (
            solutionVariablesMemory,
            1
        );

    activeTensionModel->setCouplingSignalProvider(*ionicModel);
    activeTensionModel->validateProvider();


    word ionicModelName;
    solutionVariablesMemory.lookup("ionicModel") >> ionicModelName;


    fileName ionicModelOutFile =
        runTime.path()
      / (
            ionicModelName + "_" + ionicModel->tissueName() + "_" +
            stimulusIO::protocolSuffix(solutionVariablesMemory) + ".txt"
        );

    OFstream ionicOutput(ionicModelOutFile);
    ionicOutput.setf(std::ios::fixed);
    ionicOutput.precision(7);

    // Extract the names of the fields to be exported
    const wordList ionicExportNames = ionicModel->exportedFieldNames();
    if (!ionicExportNames.empty())
    {
        Info<< "Exporting fields: " << ionicExportNames << nl;
    }
    const wordList ionicDebugNames = ionicModel->debugPrintedNames();
    if (!ionicDebugNames.empty())
    {
        Info<< "Debug printing fields: " << ionicDebugNames << nl;
    }

    ionicModel->writeHeader(ionicOutput);



    word activeTensionModelName;
    solutionVariablesMemory.lookup("activeTensionModel") >> activeTensionModelName;


    fileName activeTensionOutFile =
        runTime.path()
      / (
           activeTensionModelName + "_activeTension.txt" 
        );

    OFstream taOutput(activeTensionOutFile);
    taOutput.setf(std::ios::fixed);
    taOutput.precision(7);
    // Extract the names of the active tension fields to be exported
    const wordList activeTensionExportNames = activeTensionModel->exportedFieldNames();
    if (!activeTensionExportNames.empty())
    {
        Info<< "Exporting fields: " << activeTensionExportNames << nl;
    }
    const wordList activeTensionDebugNames = activeTensionModel->debugPrintedNames();
    if (!activeTensionDebugNames.empty())
    {
        Info<< "Debug printing fields: " << activeTensionDebugNames << nl;
    }

    activeTensionModel->writeHeader(taOutput);
    

    // Loop for the ODE solver
    scalarField lambda(1, 1.0);
    scalarField Ta(1, 0.0);

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;

        Field<Field<scalar>> dummyStates(1);
        dummyStates[0].setSize(ionicModel->nEqns());
        scalarField dummyVmField(1, 0);
        scalarField dummyIonicCurrentField(1, 0);
        scalar t0 = runTime.value() - runTime.deltaTValue();
        scalar dt = runTime.deltaTValue();
        scalar t1 = t0 + dt;

        ionicModel->solveODE
        (
            t0,
            dt,
            dummyVmField, dummyIonicCurrentField,
            dummyStates
        );
        activeTensionModel->calculateTension
        (
            runTime.value(),
            dt,
            lambda,
            Ta
        );

        if (stimulusIO::shouldWriteStep(t0, t1, solutionVariablesMemory, false))
        {
            ionicModel->write(runTime.value(), ionicOutput);
            activeTensionModel->write(runTime.value(), taOutput);
        }
    }

    Info<< nl;
    runTime.printExecutionTime(Info);
    Info<< "End" << nl;
    Info << "Results written to: " << nl; 
    Info<< "1. Ionic Model: " << ionicOutput.name() << nl;
    Info<< "Format: [Time STATES ALGEBRAIC RATES]" << nl;
    Info<< "2. Active Tension: " << taOutput.name() << nl;
    Info<< "Format: [Time STATES ALGEBRAIC RATES]" << nl;

    return 0;

}
