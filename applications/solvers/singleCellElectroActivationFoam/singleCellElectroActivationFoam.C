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
    This solver is a driver for a single cell electrophysiology model from
    cellml.org. 
    The user provides an external stimulus through the stimulus protocol dictionary.

Author
    Philip Cardiff, UCD.
    Simao Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "ionicModel.H"
#include "ionicModelIO.H"

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
            true // solve Vm equation within ODE system
        );

    word modelName;
    solutionVariablesMemory.lookup("ionicModel") >> modelName;

    fileName outFile =
        ionicModelIO::createOutputFile
        (
            solutionVariablesMemory,modelName,
            ionicModel->tissueName(), runTime
        );

    OFstream output(outFile);
    ionicModel->writeHeader(output);



    //Loop for the ODE solver

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;

        Field<Field<scalar>> dummyStates(1);
        dummyStates[0].setSize(ionicModel->nEqns());
        scalarField dummyVmField(1, 0);
        scalarField dummyIonicCurrentField(1, 0);


        ionicModel->solveODE
        (
            runTime.value() - runTime.deltaTValue(),
            runTime.deltaTValue(),
            dummyVmField, dummyIonicCurrentField,
            dummyStates
        );

        ionicModelIO::writeTimestep
            (*ionicModel, output,solutionVariablesMemory, runTime);
    }

    Info<< nl;
    runTime.printExecutionTime(Info);
    Info<< "End" << nl;
    Info<< "Results written to: " << output.name() << nl;
    Info<< "Format: [Time STATES ALGEBRAIC RATES]" << nl;

    return 0;
}


