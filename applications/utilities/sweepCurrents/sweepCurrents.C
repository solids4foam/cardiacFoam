/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

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
#include "IFstream.H"
#include "DynamicList.H"
#include "stringOps.H"
#include "OSspecific.H"
#include "ionicModel.H"

using namespace Foam;

namespace
{

bool containsWord(const wordList& words, const word& target)
{
    forAll(words, i)
    {
        if (words[i] == target)
        {
            return true;
        }
    }

    return false;
}


wordList readSweepHeaderColumns(const fileName& sweepFile)
{
    IFstream inFile(sweepFile);

    if (!inFile.good())
    {
        FatalErrorInFunction
            << "Could not read sweep output file: " << sweepFile
            << exit(FatalError);
    }

    string header;
    inFile.getLine(header);

    if (header.empty())
    {
        FatalErrorInFunction
            << "Sweep output file has empty header: " << sweepFile
            << exit(FatalError);
    }

    forAll(header, i)
    {
        if (header[i] == ',')
        {
            header[i] = ' ';
        }
    }

    const auto split = stringOps::splitSpace(header);
    wordList columns(split.size());

    forAll(columns, i)
    {
        columns[i] = split[i];
    }

    return columns;
}

} // End unnamed namespace

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

    DynamicList<word> invalidCurrents;
    forAll(currents, i)
    {
        if (!containsWord(availableCurrents, currents[i]))
        {
            invalidCurrents.append(currents[i]);
        }
    }

    if (!invalidCurrents.empty())
    {
        wordList invalidList(invalidCurrents.size());
        forAll(invalidList, i)
        {
            invalidList[i] = invalidCurrents[i];
        }

        FatalErrorInFunction
            << "Requested current(s) not available for model " << modelName
            << ": " << invalidList << nl
            << "Available currents: " << availableCurrents << nl
            << "Use listCurrentsOnly true; to inspect valid options."
            << exit(FatalError);
    }

    scalar Vmin = readScalar(sweepDict.lookup("Vmin"));
    scalar Vmax = readScalar(sweepDict.lookup("Vmax"));
    label  nPts = readLabel(sweepDict.lookup("points"));
    word outputExtension =
        sweepDict.lookupOrDefault<word>("outputExtension", "txt");
    const bool printCurrentVariables =
        bool(sweepDict.lookupOrDefault<Switch>("printCurrentVariables", true));

    if (!outputExtension.empty() && outputExtension[0] == '.')
    {
        outputExtension = outputExtension.substr(1);
    }

    const fileName outputDir(runTime.path()/"postProcessing");
    mkDir(outputDir);

    for (const word& I : currents)
    {
        const fileName out
        (
            outputDir
          / (modelName + "_" + tissue + "_" + I + "_sweep." + outputExtension)
        );

        Info<< "Sweeping " << I << " → " << out << nl;

        model.sweepCurrent(I, Vmin, Vmax, nPts, out);

        if (printCurrentVariables)
        {
            const wordList columns = readSweepHeaderColumns(out);

            if (columns.empty())
            {
                Info<< "  variables for " << I << ": (none)" << nl;
                continue;
            }

            wordList deps(max(label(0), columns.size() - 1));
            for (label c = 1; c < columns.size(); ++c)
            {
                deps[c - 1] = columns[c];
            }

            Info<< "  variables for " << I << " (" << deps.size() << "): "
                << deps << nl;
        }
    }

    Info<< "\nCompleted.\n";
    return 0;
}
