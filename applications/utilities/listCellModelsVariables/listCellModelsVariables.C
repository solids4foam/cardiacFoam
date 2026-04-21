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
    listCellModelsVariables

Description
    Reads physicsProperties to determine the active cardiac physics setup and
    lists the available ionic-model variables. If an active-tension model is
    configured, its variables are listed as well. The same report is also
    written to a file in postProcessing/.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "ionicModel.H"
#include "activeTensionModel.H"
#include "OFstream.H"
#include "OSspecific.H"

using namespace Foam;

namespace
{

dictionary electroModelDict(const IOdictionary& electroDict)
{
    if (electroDict.found("electroModel"))
    {
        word electroModelName;
        electroDict.lookup("electroModel") >> electroModelName;

        const word coeffsName(electroModelName + "Coeffs");

        if (!electroDict.found(coeffsName))
        {
            FatalErrorInFunction
                << "Expected sub-dictionary '" << coeffsName
                << "' in electroProperties for electroModel '"
                << electroModelName << "'." << nl
                << "For example:" << nl
                << "  electroModel " << electroModelName << ";" << nl
                << "  " << coeffsName << nl
                << "  {" << nl
                << "      ionicModel TNNP;" << nl
                << "      tissue endocardialCells;" << nl
                << "      ..." << nl
                << "  }" << nl
                << exit(FatalError);
        }

        return electroDict.subDict(coeffsName);
    }

    // Backward-compatible path: allow a flat dictionary if no electroModel
    // key is present.
    return dictionary(electroDict);
}


dictionary activeTensionDict(const dictionary& modelDict)
{
    const dictionary& atSubDict = modelDict.subDict("activeTensionModel");
    const word atModelType(atSubDict.lookup("activeTensionModel"));

    dictionary atDict(atSubDict);
    atDict.merge(modelDict);

    // merge(modelDict) injects the parent activeTensionModel sub-dictionary,
    // which would otherwise replace this key with a dictionary entry.
    atDict.add("activeTensionModel", atModelType, true);

    return atDict;
}


void writeHeader
(
    Ostream& os,
    const word& physicsType,
    const IOdictionary& electroDict,
    const dictionary& modelDict
)
{
    os  << "\n========== listCellModelsVariables ==========\n\n"
        << "Selected physicsModel: " << physicsType << nl;

    if (electroDict.found("electroModel"))
    {
        os << "Selected electroModel: "
           << word(electroDict.lookup("electroModel")) << nl;
    }
    else
    {
        os << "Selected electroModel: flat electroProperties dictionary" << nl;
    }

    os << "Selected ionicModel: "
       << word(modelDict.lookup("ionicModel")) << nl;
}


void writeIonicReport(Ostream& os, const ionicModel& model)
{
    const wordList constantsNames = model.constantVariableNames();
    const wordList states = model.stateVariableNames();
    const wordList algebraic = model.algebraicVariableNames();
    const scalarField constants = model.constantVariableValues();
    const scalarField initialStates = model.stateVariableValues(0);

    os << nl
       << "Ionic constants (" << constants.size() << ") --> initial value"
       << nl;
    forAll(constants, i)
    {
        os << "  constants [" << i << "]";
        if (i < constantsNames.size())
        {
            os << " " << constantsNames[i];
        }
        os << " --> " << constants[i] << nl;
    }

    os << nl
       << "Ionic states (" << states.size() << ") --> initial value"
       << nl;
    forAll(states, i)
    {
        const scalar initialValue =
            (i < initialStates.size()) ? initialStates[i] : scalar(0);

        os << "  states [" << i << "] " << states[i]
           << " --> " << initialValue
           << nl;
    }

    os << nl << "Ionic algebraic (" << algebraic.size() << ")" << nl;
    forAll(algebraic, i)
    {
        os << "  algebraic [" << i << "] " << algebraic[i] << nl;
    }
}


void writeActiveTensionReport(Ostream& os, activeTensionModel& model)
{
    const wordList constantsNames = model.constantVariableNames();
    const wordList states = model.stateVariableNames();
    const wordList algebraic = model.algebraicVariableNames();
    const wordList rates = model.rateVariableNames();
    const scalarField constants = model.constantVariableValues();
    const scalarField initialStates = model.stateVariableValues(0);
    const scalarField initialRates = model.rateVariableValues(0);

    os << nl
       << "Active tension constants (" << constants.size() << ") --> initial value"
       << nl;
    forAll(constants, i)
    {
        os << "  constants [" << i << "]";
        if (i < constantsNames.size())
        {
            os << " " << constantsNames[i];
        }
        os << " --> " << constants[i] << nl;
    }

    os << nl
       << "Active tension states (" << states.size() << ") --> initial value"
       << nl;
    forAll(states, i)
    {
        const scalar initialValue =
            (i < initialStates.size()) ? initialStates[i] : scalar(0);

        os << "  states [" << i << "] " << states[i]
           << " --> " << initialValue
           << nl;
    }

    os << nl << "Active tension algebraic (" << algebraic.size() << ")" << nl;
    forAll(algebraic, i)
    {
        os << "  algebraic [" << i << "] " << algebraic[i] << nl;
    }

    os << nl
       << "Active tension rates (" << rates.size() << ") --> initial value"
       << nl;
    forAll(rates, i)
    {
        const scalar initialValue =
            (i < initialRates.size()) ? initialRates[i] : scalar(0);

        os << "  rates [" << i << "] " << rates[i]
           << " --> " << initialValue
           << nl;
    }
}

} // End unnamed namespace

int main(int argc, char* argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary physicsDict
    (
        IOobject
        (
            "physicsProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word physicsType(physicsDict.lookup("type"));

    const bool supportsElectro =
        physicsType == "electroModel"
     || physicsType == "electroMechanicalModel";

    if (!supportsElectro)
    {
        FatalErrorInFunction
            << "Unsupported physicsProperties type '" << physicsType << "'." << nl
            << "This utility currently supports: electroModel, "
            << "and electroMechanicalModel."
            << exit(FatalError);
    }

    IOdictionary electroDict
    (
        IOobject
        (
            "electroProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dictionary modelDict(electroModelDict(electroDict));
    modelDict.add("utilities", true);

    autoPtr<ionicModel> modelPtr =
        ionicModel::New(modelDict, 1, 0.01, false);

    const ionicModel& model = modelPtr();
    const bool hasActiveTension = modelDict.found("activeTensionModel");

    const fileName outputDir(runTime.path() / "postProcessing");
    mkDir(outputDir);

    const fileName reportFile(outputDir / "listCellModelsVariables.txt");
    OFstream report(reportFile);

    writeHeader(Info, physicsType, electroDict, modelDict);
    writeHeader(report, physicsType, electroDict, modelDict);
    writeIonicReport(Info, model);
    writeIonicReport(report, model);

    if (hasActiveTension)
    {
        dictionary atDict(activeTensionDict(modelDict));

        Info<< "Selected activeTensionModel: "
            << word(atDict.lookup("activeTensionModel")) << nl;
        report << "Selected activeTensionModel: "
               << word(atDict.lookup("activeTensionModel")) << nl;

        autoPtr<activeTensionModel> activeTensionPtr =
            activeTensionModel::New(atDict, 1);

        writeActiveTensionReport(Info, activeTensionPtr());
        writeActiveTensionReport(report, activeTensionPtr());
    }
    else
    {
        Info<< nl << "No activeTensionModel configured in electroProperties."
            << nl;
        report << nl
               << "No activeTensionModel configured in electroProperties."
               << nl;
    }

    Info<< nl << "Report written to " << reportFile << nl;
    Info<< "\nCompleted.\n";
    return 0;
}


// ************************************************************************* //
