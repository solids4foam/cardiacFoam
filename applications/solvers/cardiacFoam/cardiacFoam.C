/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Application
    cardiacFoam

Description
    General solver where the solved mathematical model  is chosen at run-time.

    cardiacFoam is the same as solids4Foam, except additional libraries (e.g.
    electroModels) are linked at compile time. solids4Foam can mimic the
    behaviour of cardiacFoam if all required libraries are loaded in the
    controlDict.

Author
    Philip Cardiff, UCD.
    Simao Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "physicsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    // Create the general physics class
    autoPtr<physicsModel> physics = physicsModel::New(runTime);

    while (runTime.run())
    {
        // Update deltaT, if desired, before moving to the next step
        physics().setDeltaT(runTime);

        runTime++;

        if (physics().printInfo())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        // Solve the mathematical model
        physics().evolve();

        // Let the physics model know the end of the time-step has been reached
        physics().updateTotalFields();

        if (runTime.outputTime())
        {
            physics().writeFields(runTime);
        }

        if (physics().printInfo())
        {
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
    }

    physics().end();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
