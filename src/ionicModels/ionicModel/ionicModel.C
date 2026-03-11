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
#include "ionicVariableCompatibility.H"




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
    // Required schema:
    // outputVariables
    // {
    //   ionic
    //   {
    //     export (...);
    //     debug  (...);
    //   }
    // }
    if (dict_.found("outputVariables"))
    {
        const dictionary& outDict = dict_.subDict("outputVariables");
        if (outDict.found("ionic"))
        {
            const dictionary& ionicOut = outDict.subDict("ionic");
            if (ionicOut.found("export"))
            {
                ionicOut.lookup("export") >> variableExport_;
            }

            if (ionicOut.found("debug"))
            {
                ionicOut.lookup("debug") >> debugVarNames_;
            }
        }
    }
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

bool Foam::ionicModel::startupLogEnabled() const
{
    return dict_.lookupOrDefault<Switch>("startupLog", true);
}

Foam::label Foam::ionicModel::sampleIntegrationPoint
(
    const label nPoints
) const
{
    if (nPoints <= 0)
    {
        return 0;
    }

    const label requested = dict_.lookupOrDefault<label>("initSampleCell", 0);

    return max(label(0), min(requested, nPoints - 1));
}

void Foam::ionicModel::logExportedFieldSelection() const
{
    if (!startupLogEnabled())
    {
        return;
    }

    const wordList exportNames = exportedFieldNames();
    if (!exportNames.empty())
    {
        Info<< "Exporting fields: " << exportNames << nl;
    }
}

bool Foam::ionicModel::hasSignal(const CouplingSignal s) const
{
    const label vmIdx =
        ionicVariableCompatibility::findVmStateIndex
        (
            ioStateNames(),
            ioNumStates()
        );
    const label caiIdx =
        ionicVariableCompatibility::findCaiStateIndex
        (
            ioStateNames(),
            ioNumStates()
        );

    switch (s)
    {
        case CouplingSignal::Vm:
        {
            if (ioVmTransform())
            {
                return true;
            }

            return vmIdx >= 0;
        }
        case CouplingSignal::Cai:
        {
            return caiIdx >= 0;
        }
        default:
            return false;
    }
}

Foam::scalar Foam::ionicModel::signal
(
    const label i,
    const CouplingSignal s
) const
{
    const auto* statesPtr = ioStatesPtr();
    if (!statesPtr || statesPtr->empty())
    {
        FatalErrorInFunction
            << "Requested coupling signal "
            << static_cast<int>(s)
            << " from ionicModel, but state storage is not available."
            << abort(FatalError);
    }

    if (i < 0 || i >= statesPtr->size())
    {
        FatalErrorInFunction
            << "Requested coupling signal index i=" << i
            << " but valid integration-point range is [0, "
            << (statesPtr->size() - 1) << "]."
            << abort(FatalError);
    }

    const scalarField& state = (*statesPtr)[i];

    if (s == CouplingSignal::Vm)
    {
        const ionicModelIO::VmTransform transformVm = ioVmTransform();
        if (transformVm)
        {
            return transformVm(state);
        }

        const label vmIdx =
            ionicVariableCompatibility::findVmStateIndex
            (
                ioStateNames(),
                ioNumStates()
            );
        if (vmIdx >= 0 && vmIdx < state.size())
        {
            return state[vmIdx];
        }
    }
    else if (s == CouplingSignal::Cai)
    {
        const label caiIdx =
            ionicVariableCompatibility::findCaiStateIndex
            (
                ioStateNames(),
                ioNumStates()
            );
        if (caiIdx >= 0 && caiIdx < state.size())
        {
            return state[caiIdx];
        }
    }

    FatalErrorInFunction
        << "Requested coupling signal "
        << static_cast<int>(s)
        << " from ionicModel, but this ionic model does not provide it."
        << abort(FatalError);

    return 0.0;
}
