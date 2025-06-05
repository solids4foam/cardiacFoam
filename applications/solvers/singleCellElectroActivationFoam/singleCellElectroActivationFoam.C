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

#include <math.h>
#include "tentusscher_noble_noble_panfilov_2004.H"
#include "fvCFD.H"
#include "interpolationTable.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"

    Info<< nl << "Creating STATES array" << endl;
    scalarList STATES(17, 0.0);

    Info<< nl << "Creating CONSTANTS array" << endl;
    scalarList CONSTANTS(46, 0.0);

    Info<< nl << "Creating ALGEBRAIC array" << endl;
    scalarList ALGEBRAIC(69, 0.0);

    Info<< nl << "Creating RATES array" << endl;
    scalarList RATES(17, 0.0);

    Info<< nl << "Calling initConsts" << endl;
    initConsts(CONSTANTS.data(), RATES.data(), STATES.data());

    Info<< nl
        << "CONSTANTS = " << CONSTANTS << nl
        << "STATES = " << STATES << endl;

    // Create the voltage vs time interpolation table
    interpolationTable<scalar> timeVsVoltage("timeVsVoltage.txt");

    // Create a file for the output
    OFstream output("output.txt");

    // Loop through time
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Time in milliseconds
        const scalar VOI = runTime.value()*1000;

        // Lookup the voltage for this time (time should be in seconds!)
        const scalar voltage = timeVsVoltage(runTime.value());

        // Set the voltage in the STATES vector
        STATES[0] = voltage;

        // TODO: we need to solve the ODEs here!

        // Compute the rates
        Info<< "Calling computeRates" << endl;
        computeRates
        (
            VOI,
            CONSTANTS.data(),
            RATES.data(),
            STATES.data(),
            ALGEBRAIC.data()
        );

        // Compute the variables
        Info<< "Calling computeVariables" << endl;
        computeVariables
        (
            VOI,
            CONSTANTS.data(),
            RATES.data(),
            STATES.data(),
            ALGEBRAIC.data()
        );

        // Write the results
        Info<< "Writing the output" << endl;
        output
            << VOI;
        forAll(STATES, i)
        {
            output
                << " " << STATES[i];
        }
        forAll(ALGEBRAIC, i)
        {
            output
                << " " << ALGEBRAIC[i];
        }
        forAll(RATES, i)
        {
            output
                << " " << RATES[i];
        }
        output
            << endl;
    }

    Info<< "The results are printed to output.txt as "
        << "[Time STATES ALGEBRAIC RATES]" << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
