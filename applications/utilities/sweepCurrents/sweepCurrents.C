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
#include "Switch.H"
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
    const wordList availableCurrents = model.availableSweepCurrents();

    const bool listCurrentsOnly =
        bool(sweepDict.lookupOrDefault<Switch>("listCurrentsOnly", false))
     || bool(sweepDict.lookupOrDefault<Switch>("listCurrents", false));

    if (listCurrentsOnly)
    {
        Info<< "Available sweep currents for model " << modelName
            << " (" << availableCurrents.size() << "):" << nl;

        if (availableCurrents.empty())
        {
            Info<< "  (none available in dependency map)" << nl;
        }
        else
        {
            forAll(availableCurrents, i)
            {
                Info<< "  - " << availableCurrents[i] << nl;
            }
        }

        Info<< "\nCompleted.\n";
        return 0;
    }

    wordList currents;
    if (sweepDict.found("currents"))
    {
        sweepDict.lookup("currents") >> currents;
    }
    else
    {
        currents = availableCurrents;
    }

    if (currents.empty())
    {
        FatalErrorInFunction
            << "No currents provided for sweepCurrents." << nl
            << "Set 'currents (...)' in constant/sweepCurrents, or set"
            << " listCurrentsOnly true to print available names."
            << exit(FatalError);
    }

    scalar Vmin = readScalar(sweepDict.lookup("Vmin"));
    scalar Vmax = readScalar(sweepDict.lookup("Vmax"));
    label  nPts = readLabel(sweepDict.lookup("points"));
    word outputExtension =
        sweepDict.lookupOrDefault<word>("outputExtension", "txt");

    if (!outputExtension.empty() && outputExtension[0] == '.')
    {
        outputExtension = outputExtension.substr(1);
    }

    for (const word& I : currents)
    {
        fileName out =
            modelName + "_" + tissue +"_" + I + "_sweep." + outputExtension;

        Info<< "Sweeping " << I << " → " << out << nl;

        model.sweepCurrent(I, Vmin, Vmax, nPts, out);
    }

    Info<< "\nCompleted.\n";
    return 0;
}
