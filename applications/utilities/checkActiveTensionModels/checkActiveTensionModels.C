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
    checkActiveTensionModels

Description
    Lists available variable names for a selected active tension model
    (and their initial values when available).

Author
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "activeTensionModel.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "\n========== checkActiveTensionModels ==========\n\n";

    IOdictionary dict
    (
        IOobject
        (
            "activeTensionProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word modelType(dict.lookup("activeTensionModel"));
    Info<< "Selected activeTensionModel: " << modelType << nl;

    autoPtr<activeTensionModel> modelPtr =
        activeTensionModel::New(dict, 1);

    activeTensionModel& model = modelPtr();
    const wordList fields = model.availableFieldNames();
    scalarField values;
    model.fieldValues(0, values);

    Info<< "\nfields (" << fields.size() << ") --> initial value" << nl;
    forAll(fields, i)
    {
        const scalar initialValue =
            (i < values.size()) ? values[i] : scalar(0);
        Info<< "  field [" << i << "] " << fields[i]
            << " --> " << initialValue
            << nl;
    }

    Info<< "\nCompleted.\n";
    return 0;
}


// ************************************************************************* //
