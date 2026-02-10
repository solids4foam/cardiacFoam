/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Application
    setFibreField

Description
    Read Ionic Currents dependencies and gating variables evolution with Voltage

Author
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "ionicModel.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "\n========== sweepCurrents ==========\n\n";

    IOdictionary sweepDict
    (
        IOobject
        (
            "sweepCurrents",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    sweepDict.add("utilities", true);
    word modelName(sweepDict.lookup("ionicModel"));
    word tissue(sweepDict.lookup("tissue"));

    autoPtr<ionicModel> modelPtr =
        ionicModel::New(sweepDict, 1, 0.01, false);

    ionicModel& model = modelPtr();

    wordList currents;
    sweepDict.lookup("currents") >> currents;

    scalar Vmin = readScalar(sweepDict.lookup("Vmin"));
    scalar Vmax = readScalar(sweepDict.lookup("Vmax"));
    label  nPts = readLabel(sweepDict.lookup("points"));

    for (const word& I : currents)
    {
        fileName out = modelName + "_" + tissue +"_" + I + "_sweep.csv";

        Info<< "Sweeping " << I << " → " << out << nl;

        model.sweepCurrent(I, Vmin, Vmax, nPts, out);
    }

    Info<< "\nCompleted.\n";
    return 0;
}

