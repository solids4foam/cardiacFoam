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

\*---------------------------------------------------------------------------*/

#include "ionicModel.H"
#include "ionicSelector.H"




// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ionicModel, 0);
    defineRunTimeSelectionTable(ionicModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionicModel::ionicModel
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ODESystem(),
    odeSolver_(),
    dict_(dict),
    step_(num, initialDeltaT),
    tissue_(-1),
    solveVmWithinODESolver_(solveVmWithinODESolver)
{
    if (dict_.found("IonicModelVariables_export"))
        dict_.lookup("IonicModelVariables_export") >> variableExport_;
    else if (dict_.found("exportedVariables"))
        dict_.lookup("exportedVariables") >> variableExport_;

    if (dict_.found("IonicModelVariables_debugPrint"))
        dict_.lookup("IonicModelVariables_debugPrint") >> debugVarNames_;
    else if (dict_.found("debugPrintVariables"))
        dict_.lookup("debugPrintVariables") >> debugVarNames_;
}

void::Foam::ionicModel::setTissueFromDict()
{
    tissue_ =
        ionicSelector::selectTissueOrDimension
        (
            dict_, hasManufacturedSolution(),
            supportedTissueTypes(), supportedDimensions()
        );
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ionicModel>
Foam::ionicModel::New
(
    const dictionary& dict,
    const label nIntegrationPoints,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
{
    const word modelType(dict.lookup("ionicModel"));
    Info<< "Selecting ionic model: " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "ionicModel",
            modelType,
            *dictionaryConstructorTablePtr_
        )   << exit(FatalIOError);
    }

    return autoPtr<ionicModel>
    (
        ctorPtr(dict, nIntegrationPoints, initialDeltaT, solveVmWithinODESolver)
    );

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ionicModel::~ionicModel()
{}


// ************************************************************************* //

bool Foam::ionicModel::utilitiesMode() const
{
return dict_.found("utilities") && readBool(dict_.lookup("utilities"));
}


bool Foam::ionicModel::hasSignal(const CouplingSignal) const
{
    // Base class provides no coupling signals by default
    return false;
}

Foam::scalar Foam::ionicModel::signal
(
    const label,
    const CouplingSignal s
) const
{
    FatalErrorInFunction
        << "Requested coupling signal "
        << static_cast<int>(s)
        << " from ionicModel, but this ionic model does not provide it."
        << nl
        << "You must override hasSignal() and signal() in the derived ionic model."
        << abort(FatalError);

    return 0.0;
}
