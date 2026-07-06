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
#include "ionicHeterogeneity.H"
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
      solveVmWithinODESolver_(solveVmWithinODESolver)
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


void Foam::ionicModel::configureApexBaseBandsHeterogeneity
(
    const scalarField& apexDist,
    const dictionary& dict
)
{
    (void)apexDist;
    (void)dict;

    FatalErrorInFunction
        << "apexBaseBands heterogeneity was requested for ionic model " << type()
        << ", but this model does not support apex-to-base heterogeneity."
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

void Foam::ionicModel::configureApexBaseBandsHeterogeneityImpl
(
    const scalarField& apexDist,
    const dictionary& dict,
    PtrList<scalarField>& heterogeneousConstants
) const
{
    const scalar beta =
        dict.lookupOrDefault<scalar>("beta", 3.0);
    const scalar scalingMin =
        dict.lookupOrDefault<scalar>("scalingMin", 0.2);
    const scalar scalingMax =
        dict.lookupOrDefault<scalar>("scalingMax", 5.0);
    const wordList variables(dict.lookup("variables"));

    if (variables.empty())
    {
        FatalErrorInFunction
            << "apexBaseBands: 'variables' list is empty for ionic model "
            << type() << ". Specify at least one constant name to scale."
            << exit(FatalError);
    }

    if (scalingMin <= 0.0 || scalingMax <= 0.0 || scalingMax < scalingMin)
    {
        FatalErrorInFunction
            << "apexBaseBands: invalid scalingMin=" << scalingMin
            << " scalingMax=" << scalingMax
            << ". Require 0 < scalingMin <= scalingMax."
            << exit(FatalError);
    }

    const label nConst = ioNumConstants();
    const char* const* names = ioConstantNames();
    const scalarField* baseConstants = ioConstantsPtr();

    if (!names || nConst <= 0 || !baseConstants || baseConstants->empty())
    {
        FatalErrorInFunction
            << "apexBaseBands was requested for ionic model " << type()
            << ", but this model does not expose constant metadata "
            << "(ioConstantNames / ioConstantsPtr)."
            << exit(FatalError);
    }

    labelList indices(variables.size(), -1);
    forAll(variables, vi)
    {
        for (label ci = 0; ci < nConst; ci++)
        {
            if (word(names[ci]) == variables[vi])
            {
                indices[vi] = ci;
                break;
            }
        }
        if (indices[vi] < 0)
        {
            FatalErrorInFunction
                << "apexBaseBands: variable '" << variables[vi]
                << "' not found in constant names of ionic model " << type()
                << ". Available constants: ";
            for (label ci = 0; ci < nConst; ci++)
            {
                FatalErrorInFunction << names[ci] << ' ';
            }
            FatalErrorInFunction << exit(FatalError);
        }
    }

    if (heterogeneousConstants.empty())
    {
        heterogeneousConstants.setSize(apexDist.size());
        forAll(apexDist, cellI)
        {
            heterogeneousConstants.set(cellI, new scalarField(*baseConstants));
        }
    }
    else if (heterogeneousConstants.size() != apexDist.size())
    {
        FatalErrorInFunction
            << "apexBaseBands: longitudinal distance field has "
            << apexDist.size() << " values but " << type()
            << " has " << heterogeneousConstants.size()
            << " heterogeneous constant sets."
            << exit(FatalError);
    }

    forAll(apexDist, cellI)
    {
        const scalar d = min(max(apexDist[cellI], scalar(0.0)), scalar(1.0));
        const scalar f =
            ionicHeterogeneity::apexBaseScale(d, beta, scalingMin, scalingMax);

        scalarField& consts = heterogeneousConstants[cellI];
        forAll(indices, vi)
        {
            consts[indices[vi]] *= f;
        }
    }
}


void Foam::ionicModel::configureTransmuralBandHeterogeneity
(
    const scalarField& transmuralDistance,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
) const
{
    const auto* statesPtr = ioStatesPtr();

    if (statesPtr && transmuralDistance.size() != statesPtr->size())
    {
        FatalErrorInFunction
            << "Transmural distance field has " << transmuralDistance.size()
            << " values, but " << type() << " was configured with "
            << statesPtr->size() << " integration points."
            << exit(FatalError);
    }

    const word mode =
        heterogeneityDict.lookupOrDefault<word>("mode", "transmuralBands");

    if (mode == "namedRegions")
    {
        configureNamedRegionHeterogeneity
        (
            transmuralDistance, heterogeneityDict, heterogeneousConstants,
            heterogeneousInitialStates
        );
        return;
    }

    if (mode != "transmuralBands")
    {
        FatalErrorInFunction
            << "Unsupported " << type() << " ionicHeterogeneity mode '"
            << mode << "'. Supported modes: transmuralBands, namedRegions."
            << exit(FatalError);
    }

    const word smoothing =
        heterogeneityDict.lookupOrDefault<word>("smoothing", "smoothstep");
    const word transitionMode =
        heterogeneityDict.lookupOrDefault<word>("transitionMode", "blend");
    const scalar endoMInterface =
        heterogeneityDict.lookupOrDefault<scalar>("endoMInterface", 0.3);
    const scalar mEpiInterface =
        heterogeneityDict.lookupOrDefault<scalar>("mEpiInterface", 0.7);
    const scalar transitionWidth =
        heterogeneityDict.lookupOrDefault<scalar>("transitionWidth", 0.1);

    ionicHeterogeneity::validateTransmuralBandConfig
    (
        endoMInterface,
        mEpiInterface,
        transitionWidth,
        smoothing,
        transitionMode
    );

    const scalarField endoConstants =
        constantsForTissue(ionicSelector::tissueFlag("endocardialCells"));
    const scalarField mCellConstants =
        constantsForTissue(ionicSelector::tissueFlag("mCells"));
    const scalarField epiConstants =
        constantsForTissue(ionicSelector::tissueFlag("epicardialCells"));

    if
    (
        endoConstants.empty()
     || endoConstants.size() != mCellConstants.size()
     || endoConstants.size() != epiConstants.size()
    )
    {
        FatalErrorInFunction
            << "ionicHeterogeneity was requested for ionic model " << type()
            << ", but this model does not provide endo/M/epi constants."
            << exit(FatalError);
    }

    scalarField endoStates, mCellStates, epiStates;
    bool blendStates = false;
    if (heterogeneousInitialStates)
    {
        endoStates = initialStatesForTissue(ionicSelector::tissueFlag("endocardialCells"));
        mCellStates = initialStatesForTissue(ionicSelector::tissueFlag("mCells"));
        epiStates = initialStatesForTissue(ionicSelector::tissueFlag("epicardialCells"));

        if
        (
            !endoStates.empty()
         && endoStates.size() == mCellStates.size()
         && endoStates.size() == epiStates.size()
        )
        {
            blendStates = true;
            heterogeneousInitialStates->clear();
            heterogeneousInitialStates->setSize(transmuralDistance.size());
        }
    }

    heterogeneousConstants.clear();
    heterogeneousConstants.setSize(transmuralDistance.size());

    forAll(transmuralDistance, integrationPtI)
    {
        const scalar rawT = transmuralDistance[integrationPtI];

        if (rawT < -SMALL || rawT > 1.0 + SMALL)
        {
            FatalErrorInFunction
                << "Transmural distance value t=" << rawT
                << " at integration point " << integrationPtI
                << " is outside the expected [0, 1] range."
                << exit(FatalError);
        }

        const scalar t = min(max(rawT, scalar(0.0)), scalar(1.0));
        scalarField mappedConstants(endoConstants.size(), 0.0);
        const ionicHeterogeneity::TransmuralBandWeights weights =
            ionicHeterogeneity::transmuralBandWeights
            (
                t,
                endoMInterface,
                mEpiInterface,
                transitionWidth,
                smoothing,
                transitionMode
            );

        forAll(mappedConstants, constantI)
        {
            mappedConstants[constantI] =
                weights.endo*endoConstants[constantI]
              + weights.mCell*mCellConstants[constantI]
              + weights.epi*epiConstants[constantI];
        }

        heterogeneousConstants.set
        (
            integrationPtI,
            new scalarField(mappedConstants)
        );

        if (blendStates)
        {
            scalarField mappedStates(endoStates.size(), 0.0);
            forAll(mappedStates, stateI)
            {
                mappedStates[stateI] =
                    weights.endo*endoStates[stateI]
                  + weights.mCell*mCellStates[stateI]
                  + weights.epi*epiStates[stateI];
            }
            heterogeneousInitialStates->set
            (
                integrationPtI,
                new scalarField(mappedStates)
            );
        }
    }


}


void Foam::ionicModel::configureNamedRegionHeterogeneity
(
    const scalarField& fieldValues,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
) const
{
    const word smoothing =
        heterogeneityDict.lookupOrDefault<word>("smoothing", "smoothstep");
    const word transitionMode =
        heterogeneityDict.lookupOrDefault<word>("transitionMode", "blend");
    const scalar transitionWidth =
        heterogeneityDict.lookupOrDefault<scalar>("transitionWidth", 0.1);

    if (transitionMode != "blend" && transitionMode != "hard")
    {
        FatalErrorInFunction
            << "Unsupported ionicHeterogeneity transitionMode '"
            << transitionMode << "' for mode namedRegions. Supported: "
            << "blend, hard."
            << exit(FatalError);
    }

    if (!heterogeneityDict.found("regions"))
    {
        FatalErrorInFunction
            << "ionicHeterogeneity mode namedRegions requires a 'regions' "
            << "sub-dictionary for ionic model " << type() << "."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedFieldRegion> regions =
        ionicHeterogeneity::parseNamedFieldRegions
        (
            heterogeneityDict.subDict("regions")
        );

    wordList regionNames(regions.size());
    forAll(regions, i)
    {
        regionNames[i] = regions[i].name;
    }

    PtrList<scalarField> regionConstants(regions.size());
    PtrList<scalarField> regionStates(regions.size());
    bool blendStates = (heterogeneousInitialStates != nullptr);

    forAll(regions, i)
    {
        regionConstants.set
        (
            i,
            new scalarField
            (
                constantsForRegion(regionNames[i], regions[i].baseline, regionNames)
            )
        );

        if (regionConstants[i].empty())
        {
            FatalErrorInFunction
                << "ionicHeterogeneity mode namedRegions was requested for "
                << "ionic model " << type() << ", but constantsForTissue() "
                << "returned no constants. This model does not support "
                << "region-based heterogeneity."
                << exit(FatalError);
        }

        if (blendStates)
        {
            regionStates.set
            (
                i,
                new scalarField
                (
                    initialStatesForRegion
                    (
                        regionNames[i], regions[i].baseline, regionNames
                    )
                )
            );

            if (regionStates[i].empty())
            {
                blendStates = false;
            }
        }
    }

    heterogeneousConstants.clear();
    heterogeneousConstants.setSize(fieldValues.size());

    if (blendStates)
    {
        heterogeneousInitialStates->clear();
        heterogeneousInitialStates->setSize(fieldValues.size());
    }

    forAll(fieldValues, cellI)
    {
        const scalar rawT = fieldValues[cellI];

        if (rawT < -SMALL || rawT > 1.0 + SMALL)
        {
            FatalErrorInFunction
                << "Named-region field value t=" << rawT
                << " at integration point " << cellI
                << " is outside the expected [0, 1] range."
                << exit(FatalError);
        }

        const scalar t = min(max(rawT, scalar(0.0)), scalar(1.0));

        const List<ionicHeterogeneity::NamedRegionWeight> weights =
            ionicHeterogeneity::namedRegionWeightsAt
            (
                t, regions, transitionWidth, smoothing, transitionMode
            );

        scalarField mappedConstants(regionConstants[0].size(), 0.0);
        scalarField mappedStates;
        if (blendStates)
        {
            mappedStates.setSize(regionStates[0].size(), 0.0);
        }

        forAll(weights, wI)
        {
            const label regionI = regionNames.find(weights[wI].name);
            mappedConstants += weights[wI].weight*regionConstants[regionI];

            if (blendStates)
            {
                mappedStates += weights[wI].weight*regionStates[regionI];
            }
        }

        heterogeneousConstants.set(cellI, new scalarField(mappedConstants));

        if (blendStates)
        {
            heterogeneousInitialStates->set
            (
                cellI, new scalarField(mappedStates)
            );
        }
    }
}


Foam::scalarField Foam::ionicModel::constantsForRegion
(
    const word& regionName,
    const word& baseline,
    const wordList& knownRegionNames
) const
{
    const label tissueFlag = ionicSelector::tissueFlag(baseline);

    scalarField constants = constantsForTissue(tissueFlag);

    if (constants.empty())
    {
        return constants;
    }

    ionicModelIO::applyConstantOverrides
    (
        constants,
        ioConstantNames(),
        ioNumConstants(),
        dict_,
        type(),
        regionName,
        knownRegionNames
    );

    return constants;
}


Foam::scalarField Foam::ionicModel::initialStatesForRegion
(
    const word& regionName,
    const word& baseline,
    const wordList& knownRegionNames
) const
{
    (void)regionName;
    (void)knownRegionNames;

    const label tissueFlag = ionicSelector::tissueFlag(baseline);

    return initialStatesForTissue(tissueFlag);
}

