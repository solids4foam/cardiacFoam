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
    cellml.org. In this case, it is used to drive the Ten Tusscher, Noble,
    Noble, Panfilov, 2004 model.

    The user provides the transmembrane voltage vs time as an input.

Author
    Philip Cardiff, UCD.
    Simao Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ionicModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createFields.H"


    // Create ionicModelCellML object
    autoPtr<ionicModel> ionicModel =
    ionicModel::New
    (
        solutionVariablesMemory,
        1, //one integration point
        runTime.deltaTValue(),
        true // solve Vm equation within ODE system
    );
    word modelName;
    solutionVariablesMemory.lookup("ionicModel") >> modelName;
    word tissueName;
    solutionVariablesMemory.lookup("tissue") >> tissueName;
    // Read pacing parameters
    label s1, s2;
    solutionVariablesMemory.lookup("stim_period_S1") >> s1;
    solutionVariablesMemory.lookup("stim_period_S2") >> s2;
    // Detect protocol
    bool protocolMode = (s2 > 0);
    // Build output filename
    Foam::fileName modelOutputFile;
    if (protocolMode)
    {
        // S1-S2 protocol naming
        Foam::word s1Str = Foam::name(s1);
        Foam::word s2Str = Foam::name(s2);

        modelOutputFile =
            runTime.path()
            / (modelName + "_"
            + tissueName + "_"
            + s1Str + "-" + s2Str + "ms.txt");

        Info << "Detected S1-S2 protocol: S1=" << s1
            << " ms, S2=" << s2 << " ms" << endl;
    }
    else
    {
        // Regular pacing naming
        modelOutputFile =
            runTime.path()
            / (modelName + "_" + tissueName+ ".txt");

        Info << "Detected periodic pacing: BCL=" << s1 << " ms" << endl;
    }

    // Open output file
    Foam::OFstream output(modelOutputFile);

    // Write header
    ionicModel->writeHeader(output);

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;

        // Allocate states storage for 1 cell
        Field<Field<scalar>> dummyStates(1);
        dummyStates[0].setSize(ionicModel->nEqns());
        scalarField dummyVmField(1, 0);
        scalarField dummyIonicCurrentField(1, 0);

        // Solve the ionic model given the current voltage and calculate the
        // ionic model currents
        ionicModel->solveODE
        (
            runTime.value() - runTime.deltaTValue(),
            runTime.deltaTValue(),
            dummyVmField,
            dummyIonicCurrentField,
            dummyStates
        );

        if (protocolMode)
        {
            // S1–S2 protocol normally writes after a long warm-up window
            if (runTime.value() > 8)
            {
                ionicModel->write(runTime.value(), output);
            }
        }
        else
        {
            ionicModel->write(runTime.value(), output);
        }
    }
    
    Info<< nl << endl;
    runTime.printExecutionTime(Info);
    Info<< "End" << nl << endl;
    Info<< "Results written to: " << output.name() << nl;
    Info<< "Format: [Time STATES ALGEBRAIC RATES]" << nl;

    return 0;
}


// ************************************************************************* //
