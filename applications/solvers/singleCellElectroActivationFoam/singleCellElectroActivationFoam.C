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

\*---------------------------------------------------------------------------*/

// #include <math.h>
// #include "tentusscher_noble_noble_panfilov_2004.H"
#include "fvCFD.H"
#include "interpolationTable.H"
#include "OFstream.H"
#include <string>
#include "ionicModelCellML.H"        // Base class
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"

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

    // Model selection (e.g., Gaur, TNNP, Courtemanche)
    
    word modelName;
    electroActivationProperties.lookup("ionicModel") >> modelName;
    
    Info<< "Using ionic model: " << modelName << nl << endl;

    word tissueName;
    electroActivationProperties.lookup("tissue") >> tissueName;

    std::string fileName = modelName + "_" + tissueName + "_output.txt";
    OFstream output(fileName);

   

    // Create ionicModelCellML object
    ionicModelCellML ionicModel(electroActivationProperties);
    
    // Write header once
    ionicModel.writeHeader(output);

    // Loop through time
    Info<< "\nStarting time loop\n" << endl;
    
        while (runTime.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
    
            // Time in milliseconds
            const scalar VOI = runTime.value()*1000;
    
            // Define the time step in ms
            const scalar deltaT = runTime.deltaTValue()*1000;
    
            // Define the old time in ms
            const scalar tOld = VOI - deltaT;
    
            // Solve the ionic model ODEs for this time step given the known voltage
            ionicModel.solve(tOld, deltaT);
    
            // Compute the variables for the given time
            ionicModel.computeVariables(VOI);
    
            // Write the fields to a file
            ionicModel.write(VOI, output);
        }

    Info<< "Simulation complete.\n"
        << "Results written to: " << output.name() << nl
        << "Format: [Time STATES ALGEBRAIC RATES]" << nl;

    return 0;
}

// ************************************************************************* //