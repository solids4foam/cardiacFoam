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

#include "myocardiumPrePacing.H"
#include "ionicHeterogeneity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myocardiumPrePacing, 0);
}

namespace
{

//- Keep the entry of the processor that computed it
struct takeComputed
{
    void operator()(Foam::scalarField& x, const Foam::scalarField& y) const
    {
        if (x.empty())
        {
            x = y;
        }
    }
};

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myocardiumPrePacing::myocardiumPrePacing
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    const scalar dt,
    const scalarField* heterogeneityFieldValues,
    const label nCells
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    electroProperties_(electroProperties),
    dt_(dt),
    regionNames_
    (
        1, electroProperties.lookupOrDefault<word>("tissue", "myocyte")
    ),
    regionFieldValue_(1, 0.0),
    cellRegion_(nCells, 0),
    nBlendedCells_(0),
    pureHeterogeneityDict_(),
    hasHeterogeneity_(electroProperties.found("ionicHeterogeneity")),
    configs_(),
    regionPaced_(),
    ionicStates_()
{
    if (hasHeterogeneity_)
    {
        const dictionary& heterogeneityDict =
            electroProperties.subDict("ionicHeterogeneity");
        const word mode
        (
            heterogeneityDict.lookupOrDefault<word>("mode", word::null)
        );

        if (!heterogeneityFieldValues || heterogeneityFieldValues->size() != nCells)
        {
            FatalErrorInFunction
                << "prePacing needs the ionicHeterogeneity field value of "
                << "all " << nCells << " myocardium cells."
                << exit(FatalError);
        }
        const scalarField& fieldValues = *heterogeneityFieldValues;

        pureHeterogeneityDict_ = heterogeneityDict;

        if (mode == "cellZoneRegions")
        {
            const List<ionicHeterogeneity::NamedCellZoneRegion> regions =
                ionicHeterogeneity::parseNamedCellZoneRegions
                (
                    heterogeneityDict.subDict("regions")
                );

            regionNames_.setSize(regions.size());
            regionFieldValue_.setSize(regions.size());
            forAll(regions, regionI)
            {
                regionNames_[regionI] = regions[regionI].name;
                regionFieldValue_[regionI] = regionI;
            }

            forAll(cellRegion_, cellI)
            {
                cellRegion_[cellI] = label(round(fieldValues[cellI]));
            }
        }
        else if (mode == "namedRegions")
        {
            const List<ionicHeterogeneity::NamedFieldRegion> regions =
                ionicHeterogeneity::parseNamedFieldRegions
                (
                    heterogeneityDict.subDict("regions")
                );

            regionNames_.setSize(regions.size());
            regionFieldValue_.setSize(regions.size());
            forAll(regions, regionI)
            {
                regionNames_[regionI] = regions[regionI].name;
                regionFieldValue_[regionI] =
                    0.5*(regions[regionI].rangeMin + regions[regionI].rangeMax);
            }

            const word transitionMode
            (
                heterogeneityDict.lookupOrDefault<word>
                (
                    "transitionMode", "hard"
                )
            );
            const scalar transitionWidth =
                heterogeneityDict.lookupOrDefault<scalar>
                (
                    "transitionWidth", 0.0
                );
            const word smoothing
            (
                heterogeneityDict.lookupOrDefault<word>
                (
                    "smoothing", "smoothstep"
                )
            );

            forAll(cellRegion_, cellI)
            {
                const scalar t =
                    min(max(fieldValues[cellI], scalar(0.0)), scalar(1.0));
                const List<ionicHeterogeneity::NamedRegionWeight> weights =
                    ionicHeterogeneity::namedRegionWeightsAt
                    (
                        t, regions, transitionWidth, smoothing, transitionMode
                    );

                label dominant = 0;
                forAll(weights, wI)
                {
                    if (weights[wI].weight > weights[dominant].weight)
                    {
                        dominant = wI;
                    }
                }

                cellRegion_[cellI] = regionNames_.find(weights[dominant].name);
                if (weights.size() > 1)
                {
                    ++nBlendedCells_;
                }
            }

            // The single cell sits mid-range with no transition blend.
            pureHeterogeneityDict_.set("transitionMode", word("hard"));
        }
        else
        {
            FatalErrorInFunction
                << "constant/prePacingProperties enables prePacing, but "
                << "ionicHeterogeneity mode '" << mode << "' is not "
                << "supported. Supported: namedRegions, cellZoneRegions."
                << exit(FatalError);
        }

        if (heterogeneityDict.found("gradientAxes"))
        {
            Info<< "prePacing: gradientAxes scaling is applied to the "
                << "tissue only; cells start from their main tissue's "
                << "pre-paced state." << endl;
        }
    }

    configs_.setSize(regionNames_.size());
    forAll(regionNames_, regionI)
    {
        configs_[regionI] = prePacingIO::configFor(mesh_, regionNames_[regionI]);
    }

    regionPaced_.setSize(regionNames_.size(), false);
    ionicStates_.setSize(regionNames_.size());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::myocardiumPrePacing::singleCellDt
(
    const label regionI
) const
{
    const scalar deltaT = config(regionI).deltaT;
    return deltaT > 0 ? deltaT : dt_;
}


Foam::autoPtr<Foam::ionicModel> Foam::myocardiumPrePacing::newSingleCell
(
    const label regionI,
    const bool fromPacedState
) const
{
    const prePacingIO::PrePacingConfig& cfg = config(regionI);
    const word& regionName = regionNames_[regionI];

    // Pacing protocol: prePacingProperties' singleCellStimulus, else the
    // ionic dictionary's own; nstim1 is forced to maxBeats so every beat is
    // paced (its default is a single pulse).
    dictionary modelDict(electroProperties_);
    if (cfg.hasSingleCellStimulus)
    {
        modelDict.set("singleCellStimulus", cfg.singleCellStimulus);
    }

    if (modelDict.isDict("singleCellStimulus"))
    {
        dictionary stimDict(modelDict.subDict("singleCellStimulus"));
        stimDict.set("nstim1", cfg.maxBeats);
        modelDict.set("singleCellStimulus", stimDict);
    }
    else if (!cfg.beatComparisonIntervalSet)
    {
        FatalErrorInFunction
            << "constant/prePacingProperties enables prePacing for region '"
            << regionName << "', but no pacing protocol is defined: add a "
            << "'singleCellStimulus' sub-dictionary (stim_start, "
            << "stim_period_S1, stim_duration, stim_amplitude; ms) to "
            << "constant/prePacingProperties. For a self-beating model left "
            << "unpaced, set 'beatComparisonInterval' (ms, its own beat "
            << "period) explicitly instead."
            << exit(FatalError);
    }

    if (!cfg.singleCellIonicModel.empty())
    {
        modelDict.set("ionicModel", cfg.singleCellIonicModel);
    }
    if (!cfg.batchedIntegrator.empty())
    {
        modelDict.set("batchedIntegrator", cfg.batchedIntegrator);
    }

    autoPtr<ionicModel> modelPtr =
        ionicModel::New(modelDict, 1, singleCellDt(regionI), true);

    if (hasHeterogeneity_)
    {
        modelPtr->configureIonicHeterogeneity
        (
            scalarField(1, regionFieldValue_[regionI]),
            pureHeterogeneityDict_
        );
    }

    const PtrList<scalarField>* statesPtr = modelPtr->ioStatesPtr();

    if (!statesPtr || statesPtr->size() != 1)
    {
        FatalErrorInFunction
            << "Ionic model '" << modelPtr->type() << "' does not expose a "
            << "single-cell state vector; prePacing cannot use it. Remove "
            << "constant/prePacingProperties."
            << exit(FatalError);
    }

    if (fromPacedState)
    {
        if (!regionPaced_[regionI])
        {
            FatalErrorInFunction
                << "Region '" << regionName << "' has not been pre-paced."
                << exit(FatalError);
        }
        modelPtr->setStates(List<scalarField>(1, ionicStates_[regionI]));
    }

    return modelPtr;
}


void Foam::myocardiumPrePacing::paceIonic()
{
    labelList beats(regionNames_.size(), 0);

    forAll(regionNames_, regionI)
    {
        const prePacingIO::PrePacingConfig& cfg = config(regionI);

        if (!cfg.enabled)
        {
            Info<< "prePacing: region '" << regionNames_[regionI]
                << "' disabled in constant/prePacingProperties; left at its "
                << "initial state." << endl;
            continue;
        }

        regionPaced_[regionI] = true;

        if (Pstream::myProcNo() != owner(regionI))
        {
            continue;
        }

        autoPtr<ionicModel> cellPtr = newSingleCell(regionI, false);

        beats[regionI] =
            cellPtr->prePaceToConvergence
            (
                singleCellDt(regionI),
                cfg.tolerance,
                cfg.minBeats,
                cfg.maxBeats,
                cfg.beatComparisonInterval
            );

        ionicStates_[regionI] = (*cellPtr->ioStatesPtr())[0];
    }

    Pstream::listCombineReduce(ionicStates_, takeComputed());
    Pstream::listCombineReduce(beats, maxEqOp<label>());

    forAll(regionNames_, regionI)
    {
        if (regionPaced_[regionI])
        {
            Info<< "prePacing: region '" << regionNames_[regionI]
                << "' converged after " << beats[regionI] << " beats";
            if (Pstream::parRun())
            {
                Info<< " (processor " << owner(regionI) << ")";
            }
            Info<< endl;
        }
    }
}


void Foam::myocardiumPrePacing::seedIonic
(
    ionicModel& tissueModel,
    volScalarField& Vm
) const
{
    const PtrList<scalarField>* statesPtr = tissueModel.ioStatesPtr();
    scalarField& VmValues = Vm.primitiveFieldRef();

    if
    (
        !statesPtr
     || statesPtr->size() != cellRegion_.size()
     || VmValues.size() != cellRegion_.size()
    )
    {
        FatalErrorInFunction
            << "Ionic model '" << tissueModel.type() << "' does not expose "
            << "states for all " << cellRegion_.size() << " myocardium "
            << "cells; prePacing cannot seed it. Remove "
            << "constant/prePacingProperties."
            << exit(FatalError);
    }

    // The single-cell model (possibly a batched twin) must share the
    // tissue model's state layout.
    forAll(regionNames_, regionI)
    {
        if (!regionPaced_[regionI])
        {
            continue;
        }

        const autoPtr<ionicModel> cellPtr = newSingleCell(regionI, false);
        const label nStates = tissueModel.ioNumStates();
        const char* const* tissueNames = tissueModel.ioStateNames();
        const char* const* cellNames = cellPtr->ioStateNames();

        bool sameLayout =
            nStates == cellPtr->ioNumStates()
         && ionicStates_[regionI].size() == nStates
         && tissueNames && cellNames;

        for (label i = 0; sameLayout && i < nStates; ++i)
        {
            sameLayout = (word(tissueNames[i]) == word(cellNames[i]));
        }

        if (!sameLayout)
        {
            FatalErrorInFunction
                << "prePacing single-cell model '" << cellPtr->type()
                << "' does not share the state layout of the tissue model '"
                << tissueModel.type() << "'; its states cannot seed the "
                << "tissue. Use a singleCellIonicModel with identical states "
                << "in constant/prePacingProperties."
                << exit(FatalError);
        }
        break;
    }

    const ionicModelIO::VmTransform transform = tissueModel.ioVmTransform();

    scalarField seedVmMv(regionNames_.size(), 0.0);
    forAll(regionNames_, regionI)
    {
        if (regionPaced_[regionI])
        {
            seedVmMv[regionI] =
                transform
              ? transform(ionicStates_[regionI])
              : ionicStates_[regionI][0];
        }
    }

    labelList regionCellCount(regionNames_.size(), 0);
    List<scalarField> states(*statesPtr);

    forAll(cellRegion_, cellI)
    {
        const label regionI = cellRegion_[cellI];
        if (regionPaced_[regionI])
        {
            states[cellI] = ionicStates_[regionI];
            VmValues[cellI] = seedVmMv[regionI]*1e-3;
            ++regionCellCount[regionI];
        }
    }

    if (!tissueModel.setStates(states))
    {
        FatalErrorInFunction
            << "Ionic model '" << tissueModel.type() << "' does not accept "
            << "seeded states; prePacing cannot seed it."
            << exit(FatalError);
    }

    Vm.correctBoundaryConditions();

    if (nBlendedCells_)
    {
        Info<< "prePacing: " << nBlendedCells_ << " cells in blend "
            << "transitions start from their dominant region's state."
            << endl;
    }

    forAll(regionNames_, regionI)
    {
        if (regionPaced_[regionI])
        {
            Info<< "prePacing: seeded region '" << regionNames_[regionI]
                << "' (" << regionCellCount[regionI] << " cells) from a "
                << "converged single-cell " << tissueModel.type()
                << " state, Vm = " << seedVmMv[regionI] << " mV." << endl;
        }
    }
}


// ************************************************************************* //
