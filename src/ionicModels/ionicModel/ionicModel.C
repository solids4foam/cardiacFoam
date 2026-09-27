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
#include "ionicHeterogeneityOrchestrator.H"
#include "ionicSelector.H"
#include "ionicVariableCompatibility.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(ionicModel, 0);
defineRunTimeSelectionTable(ionicModel, dictionary);
} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionicModel::ionicModel(const dictionary& dict,
                             const label num, const scalar initialDeltaT,
                             const Switch solveVmWithinODESolver)
    : ODESystem(), odeSolver_(), dict_(dict),
      step_(num, initialDeltaT), tissue_(-1), sex_(0),
      solveVmWithinODESolver_(solveVmWithinODESolver),
      VmRatePtr_(nullptr), activeVmRate_(0.0)
{
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

void ::Foam::ionicModel::setTissueFromDict()
{
    tissue_ = ionicSelector::selectTissue(dict_, supportedTissueTypes());
}

void ::Foam::ionicModel::setSexFromDict()
{
    const List<word> supported = supportedSexTypes();
    if (supported.empty())
    {
        sex_ = 0;
        return;
    }
    sex_ = ionicSelector::selectSex(dict_, supported);
}

void ::Foam::ionicModel::applyIonicConstantOverrides() const
{
    if (!dict_.found("ionicConstantOverrides"))
    {
        return;
    }

    scalarField* constantsPtr = ioMutableConstantsPtr();
    if (!constantsPtr)
    {
        FatalErrorInFunction
            << "ionicConstantOverrides was requested for ionic model "
            << type()
            << ", but this model does not expose mutable constants."
            << exit(FatalError);
    }

    ionicModelIO::applyConstantOverrides
    (
        *constantsPtr,
        ioConstantNames(),
        ioNumConstants(),
        dict_,
        type(),
        tissue_
    );
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ionicModel> Foam::ionicModel::New(
    const dictionary& dict, const label nIntegrationPoints,
    const scalar initialDeltaT, const Switch solveVmWithinODESolver)
{
    const word modelType(dict.lookup("ionicModel"));
    auto *ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup(dict, "ionicModel", modelType,
                             *dictionaryConstructorTablePtr_)
            << exit(FatalIOError);
    }

    return autoPtr<ionicModel>
    (
        ctorPtr(dict, nIntegrationPoints, initialDeltaT, solveVmWithinODESolver)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ionicModel::~ionicModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::ionicModel::utilitiesMode() const
{
    return dict_.found("utilities")
        && readBool(dict_.lookup("utilities"));
}

Foam::label
Foam::ionicModel::sampleIntegrationPoint(const label nPoints) const
{
    if (nPoints <= 0)
    {
        return 0;
    }

    const label requested =
        dict_.lookupOrDefault<label>("initSampleCell", 0);

    return max(label(0), min(requested, nPoints - 1));
}

bool Foam::ionicModel::hasSignal(const CouplingSignal s) const
{
    const label vmIdx = ionicVariableCompatibility::findVmStateIndex(
        ioStateNames(), ioNumStates());
    const label caiIdx = ionicVariableCompatibility::findCaiStateIndex(
        ioStateNames(), ioNumStates());

    switch (s)
    {
        case CouplingSignal::VM:
        {
            if (ioVmTransform())
            {
                return true;
            }

            return vmIdx >= 0;
        }
        case CouplingSignal::CAI:
        {
            return caiIdx >= 0;
        }
        default:
            return false;
    }
}

Foam::scalar Foam::ionicModel::signal(const label i,
                                      const CouplingSignal s) const
{
    const auto *statesPtr = ioStatesPtr();
    if (!statesPtr || statesPtr->empty())
    {
        FatalErrorInFunction
            << "Requested coupling signal " << static_cast<int>(s)
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

    if (s == CouplingSignal::VM)
    {
        const ionicModelIO::VmTransform transformVm = ioVmTransform();
        if (transformVm)
        {
            return transformVm(state);
        }

        const label vmIdx = ionicVariableCompatibility::findVmStateIndex(
            ioStateNames(), ioNumStates());
        if (vmIdx >= 0 && vmIdx < state.size())
        {
            return state[vmIdx];
        }
    }
    else if (s == CouplingSignal::CAI)
    {
        const label caiIdx = ionicVariableCompatibility::findCaiStateIndex(
            ioStateNames(), ioNumStates());
        if (caiIdx >= 0 && caiIdx < state.size())
        {
            return state[caiIdx];
        }
    }

    FatalErrorInFunction
        << "Requested coupling signal " << static_cast<int>(s)
        << " from ionicModel, but this ionic model does not provide it."
        << abort(FatalError);

    return 0.0;
}

void Foam::ionicModel::exportFields(const wordList& fieldNames,
                                    PtrList<volScalarField>& outFields) const
{
    if (!hasIOMetadata())
    {
        FatalErrorInFunction
            << "I/O metadata not available for ionic model " << typeName
            << ". Override exportFields() or provide io* hooks."
            << exit(FatalError);
    }

    const auto *statesPtr = ioStatesPtr();
    const auto *algebraicPtr = ioAlgebraicPtr();
    const auto *ratesPtr = ioRatesPtr();

    if (!statesPtr || !algebraicPtr || !ratesPtr)
    {
        FatalErrorInFunction
            << "I/O storage pointers are not available for ionic model "
            << typeName
            << ". Override exportFields() or provide io* hooks "
            << "(including ioRatesPtr)." << exit(FatalError);
    }

    ionicModelIO::exportStateFields
    (
        *statesPtr,
        *algebraicPtr,
        *ratesPtr,
        fieldNames,
        ioStateNames(),
        ioNumStates(),
        ioAlgebraicNames(),
        ioNumAlgebraic(),
        exportTransferSelectedPlanCache_,
        outFields
    );
}

void Foam::ionicModel::importFields(const volScalarField& Vm,
                                    const wordList& fieldNames,
                                    const PtrList<volScalarField>& inFields)
{
    PtrList<scalarField> *statesPtr =
        const_cast<PtrList<scalarField> *>(ioStatesPtr());

    if (!statesPtr || statesPtr->empty())
    {
        return;
    }

    const label nStates = ioNumStates();
    if (nStates <= 0 || !ioStateNames())
    {
        return;
    }

    ionicModelIO::importStateFields
    (
        *statesPtr,
        Vm,
        inFields,
        fieldNames,
        ioStateNames(),
        nStates,
        ioAlgebraicNames(),
        ioNumAlgebraic(),
        importTransferSelectedPlanCache_
    );
}


// * * * * * * * * * * * Heterogeneity Functions * * * * * * * * * * * * * * //

void Foam::ionicModel::configureIonicHeterogeneity
(
    const scalarField& transmuralDistance,
    const dictionary& heterogeneityDict
)
{
    (void)transmuralDistance;
    (void)heterogeneityDict;

    FatalErrorInFunction
        << "ionicHeterogeneity was requested for ionic model " << type()
        << ", but this model does not support spatial ionic heterogeneity."
        << exit(FatalError);
}


void Foam::ionicModel::configureGradientAxisHeterogeneity
(
    const scalarField& fieldValues,
    const dictionary& dict
)
{
    (void)fieldValues;
    (void)dict;

    FatalErrorInFunction
        << "gradientAxes heterogeneity was requested for ionic model " << type()
        << ", but this model does not support gradient-axis heterogeneity."
        << exit(FatalError);
}


Foam::scalarField Foam::ionicModel::constantsForTissue
(
    const label tissueFlag
) const
{
    return scalarField();
}


Foam::scalarField Foam::ionicModel::initialStatesForTissue
(
    const label tissueFlag
) const
{
    return scalarField();
}

void Foam::ionicModel::configureGradientAxisHeterogeneityImpl
(
    const scalarField& fieldValues,
    const dictionary& dict,
    PtrList<scalarField>& heterogeneousConstants
) const
{
    ionicHeterogeneityOrchestrator::configureGradientAxisHeterogeneity
    (
        *this, fieldValues, dict, heterogeneousConstants
    );
}


void Foam::ionicModel::configureRegionHeterogeneity
(
    const scalarField& transmuralDistance,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
) const
{
    ionicHeterogeneityOrchestrator::configureRegionHeterogeneity
    (
        *this, transmuralDistance, heterogeneityDict, heterogeneousConstants,
        heterogeneousInitialStates
    );
}


void Foam::ionicModel::prePaceToConvergence
(
    const scalar dt,
    const scalar tolerance,
    const label minBeats,
    const label maxBeats,
    const scalar beatComparisonInterval
)
{
    const PtrList<scalarField>* statesPtr = ioStatesPtr();

    if (!statesPtr || statesPtr->size() != 1)
    {
        FatalErrorInFunction
            << "prePaceToConvergence requires a single-integration-point "
            << "ionic model instance (nIntegrationPoints == 1) with "
            << "generic state-vector access (ioStatesPtr() != nullptr). "
            << "Got " << (statesPtr ? statesPtr->size() : -1)
            << " integration point(s) for model '" << type() << "'."
            << exit(FatalError);
    }

    // stimPeriodS1 and beatComparisonInterval are ms; t and dt are s.
    const scalar checkpoint =
        (
            (stimulusProtocol().stimPeriodS1 > SMALL)
          ? stimulusProtocol().stimPeriodS1
          : beatComparisonInterval
        )*1e-3;

    scalarField dummyVm(1, 0.0);
    scalarField dummyIm(1, 0.0);

    scalarField previous((*statesPtr)[0]);
    scalar t = 0.0;
    label consecutiveConverged = 0;

    for (label beat = 0; beat < maxBeats; ++beat)
    {
        // Absolute checkpoint times keep every comparison at the same phase.
        const scalar beatStart = t;
        const scalar beatEnd = scalar(beat + 1)*checkpoint;
        const label nFullSteps =
            label(std::floor((beatEnd - beatStart)/dt + 1e-9));

        for (label stepI = 0; stepI < nFullSteps; ++stepI)
        {
            solveODE(t, dt, dummyVm, dummyIm);
            t = beatStart + scalar(stepI + 1)*dt;
        }

        const scalar remainder = beatEnd - t;
        if (remainder > SMALL)
        {
            solveODE(t, remainder, dummyVm, dummyIm);
        }
        t = beatEnd;

        const scalarField& current = (*statesPtr)[0];

        scalar maxRelDelta = 0.0;
        forAll(current, stateI)
        {
            const scalar denom = Foam::max(mag(previous[stateI]), SMALL);
            maxRelDelta =
                Foam::max
                (
                    maxRelDelta,
                    mag(current[stateI] - previous[stateI])/denom
                );
        }

        previous = current;

        if (beat + 1 >= minBeats && maxRelDelta < tolerance)
        {
            ++consecutiveConverged;
            if (consecutiveConverged >= 2)
            {
                Info<< "prePaceToConvergence: " << type()
                    << " converged after " << (beat + 1)
                    << " beats (max relative state change "
                    << maxRelDelta << " < tolerance " << tolerance << ")"
                    << endl;
                return;
            }
        }
        else
        {
            consecutiveConverged = 0;
        }
    }

    FatalErrorInFunction
        << "prePaceToConvergence: " << type() << " did not converge within "
        << maxBeats << " beats of " << checkpoint*1e3 << " ms (tolerance "
        << tolerance << "). Increase maxBeats, loosen tolerance, or check "
        << "for sustained alternans/instability at this pacing rate in "
        << "constant/prePacingProperties."
        << exit(FatalError);
}
