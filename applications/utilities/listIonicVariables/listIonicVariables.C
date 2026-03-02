/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Application
    listIonicVariables

Description
    Lists available state and algebraic variable names for a selected
    ionic model.

Author
    Simao Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "ionicModel.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "\n========== listIonicVariables ==========\n\n";

    IOdictionary listDict
    (
        IOobject
        (
            "cardiacProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    listDict.add("utilities", true);

    autoPtr<ionicModel> modelPtr =
        ionicModel::New(listDict, 1, 0.01, false);

    const ionicModel& model = modelPtr();

    const wordList constantsNames = model.constantVariableNames();
    const wordList states = model.stateVariableNames();
    const wordList algebraic = model.algebraicVariableNames();
    const scalarField constants = model.constantVariableValues();
    const scalarField initialStates = model.stateVariableValues(0);

    Info<< "constants (" << constants.size() << ") --> initial value" << nl;
    forAll(constants, i)
    {
        Info<< "  constants [" << i << "]";
        if (i < constantsNames.size())
        {
            Info<< " " << constantsNames[i];
        }
        Info<< " --> " << constants[i] << nl;
    }

    Info<< nl;

    Info<< "states (" << states.size() << ") --> initial value" << nl;
    forAll(states, i)
    {
        const scalar initialValue =
            (i < initialStates.size()) ? initialStates[i] : scalar(0);

        Info<< "  states [" << i << "] " << states[i]
            << " --> " << initialValue
            << nl;
    }

    Info<< nl << "algebraic (" << algebraic.size() << ")" << nl;
    forAll(algebraic, i)
    {
        Info<< "  algebraic [" << i << "] " << algebraic[i] << nl;
    }

    Info<< "\nCompleted.\n";
    return 0;
}


// ************************************************************************* //
