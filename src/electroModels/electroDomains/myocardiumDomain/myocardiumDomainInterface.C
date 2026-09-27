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

#include "myocardiumDomainInterface.H"
#include "myocardiumDomain.H"
#include "eikonalMyocardiumDomain.H"
#include "ionicModel.H"
#include "ionicHeterogeneity.H"
#include "electroVerificationModel.H"
#include "error.H"
#include "fvMeshSubset.H"
#include "volFields.H"
#include "prePacingIO.H"

namespace Foam
{

namespace
{

word myocardiumSolverType(const dictionary& electroProperties)
{
    const word coeffDictName = electroProperties.dictName();

    return coeffDictName.endsWith("Coeffs")
        ? word(coeffDictName.substr(0, coeffDictName.size() - 6))
        : coeffDictName;
}


scalarField readTransmuralDistance
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    const dictionary& heterogeneityDict
)
{
    scalarField fullValues;

    const word mode(heterogeneityDict.lookup("mode"));

    if (mode == "cellZoneRegions")
    {
        fullValues.setSize(mesh.nCells(), -1.0);
        const dictionary& regionsDict = heterogeneityDict.subDict("regions");

        // Parse cellZone regions with the shared validator.
        const List<ionicHeterogeneity::NamedCellZoneRegion> regions =
            ionicHeterogeneity::parseNamedCellZoneRegions(regionsDict);

        forAll(regions, regionIndex)
        {
            const word& zoneName = regions[regionIndex].cellZone;
            const label zoneId = mesh.cellZones().findZoneID(zoneName);

            if (zoneId < 0)
            {
                FatalErrorInFunction
                    << "ionicHeterogeneity region '"
                    << regions[regionIndex].name << "' specifies cellZone '"
                    << zoneName << "' but it does not exist on mesh '"
                    << mesh.name() << "'."
                    << exit(FatalError);
            }

            const labelList& zoneCells = mesh.cellZones()[zoneId];
            forAll(zoneCells, i)
            {
                if (fullValues[zoneCells[i]] >= 0.0)
                {
                    FatalErrorInFunction
                        << "ionicHeterogeneity cellZoneRegions: cell "
                        << zoneCells[i] << " belongs to more than one "
                        << "region's cellZone ('"
                        << regions[regionIndex].name << "' and an earlier "
                        << "region both claim it)."
                        << exit(FatalError);
                }

                fullValues[zoneCells[i]] = scalar(regionIndex);
            }
        }
    }
    else
    {
        if (!heterogeneityDict.found("field"))
        {
            FatalErrorInFunction
                << "ionicHeterogeneity mode '" << mode << "' requires a "
                << "'field' entry naming the scalar field to read region "
                << "membership from."
                << exit(FatalError);
        }

        const word fieldName(heterogeneityDict.lookup("field"));

        const volScalarField transmuralField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        fullValues = transmuralField.primitiveField();
    }

    if (!electroProperties.found("cellZone"))
    {
        return fullValues;
    }

    const word cellZoneName(electroProperties.lookup("cellZone"));
    const label zoneId = mesh.cellZones().findZoneID(cellZoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << "Cannot find myocardium cellZone '" << cellZoneName
            << "' on mesh '" << mesh.name() << "'."
            << exit(FatalError);
    }

    fvMeshSubset subset(mesh);
    subset.setCellSubset(mesh.cellZones()[zoneId]);

    const labelUList& cellMap = subset.cellMap();
    scalarField mappedValues(cellMap.size(), 0.0);

    forAll(cellMap, subCellI)
    {
        mappedValues[subCellI] = fullValues[cellMap[subCellI]];
    }

    return mappedValues;
}


scalarField readNamedScalarField
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    const word& fieldName
)
{
    const volScalarField namedField
    (
        IOobject
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    scalarField fullValues = namedField.primitiveField();

    if (!electroProperties.found("cellZone"))
    {
        return fullValues;
    }

    const word cellZoneName(electroProperties.lookup("cellZone"));
    const label zoneId = mesh.cellZones().findZoneID(cellZoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << "Cannot find myocardium cellZone '" << cellZoneName
            << "' on mesh '" << mesh.name() << "'."
            << exit(FatalError);
    }

    fvMeshSubset subset(mesh);
    subset.setCellSubset(mesh.cellZones()[zoneId]);

    const labelUList& cellMap = subset.cellMap();
    scalarField mappedValues(cellMap.size(), 0.0);

    forAll(cellMap, subCellI)
    {
        mappedValues[subCellI] = fullValues[cellMap[subCellI]];
    }

    return mappedValues;
}

//- The model dictionary a region is pre-paced with: electroProperties
//  with the pacing protocol from prePacingProperties (else the region's own
//  singleCellStimulus), nstim1 forced to maxBeats so every beat is paced
//  (its default is a single pulse).
dictionary prePacingModelDict
(
    const dictionary& electroProperties,
    const prePacingIO::PrePacingConfig& cfg,
    const word& regionName
)
{
    dictionary modelDict(electroProperties);
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

    return modelDict;
}


//- Pace a single-cell copy of the model to its limit cycle and return
//  its state vector. With heterogeneityDict, the single cell is configured
//  as the tissue would be at heterogeneity field value fieldValue (region
//  index for cellZoneRegions, t for namedRegions); nullptr means none.
scalarField prePacedState
(
    const dictionary& modelDict,
    const prePacingIO::PrePacingConfig& cfg,
    const scalar dt,
    const dictionary* heterogeneityDict,
    const scalar fieldValue
)
{
    autoPtr<ionicModel> modelPtr = ionicModel::New(modelDict, 1, dt, true);

    if (heterogeneityDict)
    {
        modelPtr->configureIonicHeterogeneity
        (
            scalarField(1, fieldValue),
            *heterogeneityDict
        );
    }

    modelPtr->prePaceToConvergence
    (
        dt,
        cfg.tolerance,
        cfg.minBeats,
        cfg.maxBeats,
        cfg.beatComparisonInterval
    );

    return scalarField((*modelPtr->ioStatesPtr())[0]);
}


//- Pre-pace one single cell per main tissue (the 'tissue' entry, or each
//  ionicHeterogeneity region) and seed every cell with its main tissue's
//  ionic states and Vm.
void prePaceAndSeed
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    const scalar dt,
    ionicModel& tissueModel,
    volScalarField& Vm
)
{
    PtrList<scalarField>* statesPtr =
        const_cast<PtrList<scalarField>*>(tissueModel.ioStatesPtr());

    if (!statesPtr)
    {
        FatalErrorInFunction
            << "Ionic model '" << tissueModel.type() << "' does not "
            << "expose generic state access (ioStatesPtr()); prePacing "
            << "cannot seed it. Remove constant/prePacingProperties."
            << exit(FatalError);
    }

    scalarField& VmValues = Vm.primitiveFieldRef();

    if (statesPtr->size() != VmValues.size())
    {
        FatalErrorInFunction
            << "Ionic model '" << tissueModel.type() << "' holds "
            << statesPtr->size() << " integration points but Vm has "
            << VmValues.size() << " cells; prePacing cannot seed it."
            << exit(FatalError);
    }

    // One pre-paced single cell per main tissue: the 'tissue' entry, or
    // each ionicHeterogeneity region with its own constants and overrides.
    // Blend transitions and gradientAxes scaling are not pre-paced: those
    // cells keep their exact constants and start from their dominant
    // region's state.
    wordList regionNames
    (
        1, electroProperties.lookupOrDefault<word>("tissue", "myocyte")
    );
    scalarList regionFieldValue(1, 0.0);
    labelList cellRegion(VmValues.size(), 0);
    label nBlendedCells = 0;
    dictionary pureHeterogeneityDict;
    const dictionary* heterogeneityDictPtr = nullptr;

    if (electroProperties.found("ionicHeterogeneity"))
    {
        const dictionary& heterogeneityDict =
            electroProperties.subDict("ionicHeterogeneity");
        const word mode
        (
            heterogeneityDict.lookupOrDefault<word>("mode", word::null)
        );

        const scalarField fieldValues =
            readTransmuralDistance(mesh, electroProperties, heterogeneityDict);

        pureHeterogeneityDict = heterogeneityDict;

        if (mode == "cellZoneRegions")
        {
            const List<ionicHeterogeneity::NamedCellZoneRegion> regions =
                ionicHeterogeneity::parseNamedCellZoneRegions
                (
                    heterogeneityDict.subDict("regions")
                );

            regionNames.setSize(regions.size());
            regionFieldValue.setSize(regions.size());
            forAll(regions, regionI)
            {
                regionNames[regionI] = regions[regionI].name;
                regionFieldValue[regionI] = regionI;
            }

            forAll(cellRegion, cellI)
            {
                cellRegion[cellI] = label(round(fieldValues[cellI]));
            }
        }
        else if (mode == "namedRegions")
        {
            const List<ionicHeterogeneity::NamedFieldRegion> regions =
                ionicHeterogeneity::parseNamedFieldRegions
                (
                    heterogeneityDict.subDict("regions")
                );

            regionNames.setSize(regions.size());
            regionFieldValue.setSize(regions.size());
            forAll(regions, regionI)
            {
                regionNames[regionI] = regions[regionI].name;
                regionFieldValue[regionI] =
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

            forAll(cellRegion, cellI)
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

                cellRegion[cellI] = regionNames.find(weights[dominant].name);
                if (weights.size() > 1)
                {
                    ++nBlendedCells;
                }
            }

            // The pre-paced cell sits mid-range with no transition blend.
            pureHeterogeneityDict.set("transitionMode", word("hard"));
        }
        else
        {
            FatalErrorInFunction
                << "constant/prePacingProperties enables prePacing, but "
                << "ionicHeterogeneity mode '" << mode << "' is not "
                << "supported. Supported: namedRegions, cellZoneRegions."
                << exit(FatalError);
        }

        heterogeneityDictPtr = &pureHeterogeneityDict;

        if (heterogeneityDict.found("gradientAxes"))
        {
            Info<< "prePacing: gradientAxes scaling is applied to the "
                << "tissue only; cells start from their main tissue's "
                << "pre-paced state." << endl;
        }
    }

    const ionicModelIO::VmTransform transform = tissueModel.ioVmTransform();

    boolList regionPaced(regionNames.size(), false);
    List<scalarField> seedStates(regionNames.size());
    scalarField seedVmMv(regionNames.size(), 0.0);

    forAll(regionNames, regionI)
    {
        const word& regionName = regionNames[regionI];
        const prePacingIO::PrePacingConfig cfg =
            prePacingIO::configFor(mesh, regionName);

        if (!cfg.enabled)
        {
            Info<< "prePacing: region '" << regionName << "' disabled in "
                << "constant/prePacingProperties; left at its initial "
                << "state." << endl;
            continue;
        }

        Info<< "prePacing: pacing region '" << regionName << "'" << endl;

        seedStates[regionI] =
            prePacedState
            (
                prePacingModelDict(electroProperties, cfg, regionName),
                cfg,
                dt,
                heterogeneityDictPtr,
                regionFieldValue[regionI]
            );

        seedVmMv[regionI] =
            transform
          ? transform(seedStates[regionI])
          : seedStates[regionI][0];

        regionPaced[regionI] = true;
    }

    labelList regionCellCount(regionNames.size(), 0);

    forAll(cellRegion, cellI)
    {
        const label regionI = cellRegion[cellI];
        if (regionPaced[regionI])
        {
            (*statesPtr)[cellI] = seedStates[regionI];
            VmValues[cellI] = seedVmMv[regionI]*1e-3;
            ++regionCellCount[regionI];
        }
    }

    Vm.correctBoundaryConditions();

    if (nBlendedCells)
    {
        Info<< "prePacing: " << nBlendedCells << " cells in blend "
            << "transitions start from their dominant region's state."
            << endl;
    }

    forAll(regionNames, regionI)
    {
        if (regionPaced[regionI])
        {
            Info<< "prePacing: seeded region '" << regionNames[regionI]
                << "' (" << regionCellCount[regionI] << " cells) from a "
                << "converged single-cell " << tissueModel.type()
                << " state, Vm = " << seedVmMv[regionI] << " mV." << endl;
        }
    }
}


} // End anonymous namespace


autoPtr<myocardiumDomainInterface> myocardiumDomainInterface::New
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList& postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    autoPtr<ionicModel>& ionicModelPtr,
    autoPtr<electroVerificationModel>& verificationModelPtr,
    scalar initialDeltaT
)
{
    const word solverType = myocardiumSolverType(electroProperties);

    if (solverType == "eikonalSolver")
    {
        ionicModelPtr.clear();
        verificationModelPtr.clear();

        return autoPtr<myocardiumDomainInterface>
        (
            new eikonalMyocardiumDomain(mesh, electroProperties)
        );
    }

    ionicModelPtr =
        ionicModel::New
        (
            electroProperties,
            myocardiumDomain::configuredCellCount(mesh, electroProperties),
            initialDeltaT
        );

    if (electroProperties.found("ionicHeterogeneity"))
    {
        const dictionary& heterogeneityDict =
            electroProperties.subDict("ionicHeterogeneity");

        if (heterogeneityDict.found("field") && !heterogeneityDict.found("mode"))
        {
            FatalErrorInFunction
                << "electroProperties.ionicHeterogeneity has a 'field' "
                << "entry but no 'mode' entry. 'mode' is required whenever "
                << "region heterogeneity is configured: namedRegions or "
                << "cellZoneRegions."
                << exit(FatalError);
        }

        if (heterogeneityDict.found("mode"))
        {
            const scalarField transmuralDistance =
                readTransmuralDistance(mesh, electroProperties, heterogeneityDict);

            ionicModelPtr->configureIonicHeterogeneity
            (
                transmuralDistance,
                heterogeneityDict
            );
        }

        if (heterogeneityDict.found("gradientAxes"))
        {
            const dictionary& axesDict =
                heterogeneityDict.subDict("gradientAxes");

            forAllConstIter(dictionary, axesDict, iter)
            {
                const word axisName(iter().keyword());
                const dictionary& axisDict = axesDict.subDict(axisName);
                const word fieldName(axisDict.lookup("field"));

                const scalarField axisField =
                    readNamedScalarField(mesh, electroProperties, fieldName);

                ionicModelPtr->configureGradientAxisHeterogeneity
                (
                    axisField,
                    axisDict
                );
            }
        }
    }

    verificationModelPtr =
        electroVerificationModel::New(electroProperties);

    autoPtr<myocardiumDomain> tissuePtr
    (
        myocardiumDomain::New
        (
            mesh,
            electroProperties,
            outFields,
            postProcessFieldNames,
            postProcessFields,
            ionicModelPtr(),
            verificationModelPtr.get()
        )
    );

    // Vm exists only once myocardiumDomain::New has run -- prePacing must
    // follow it.
    if (prePacingIO::configFor(mesh, word::null).enabled)
    {
        prePaceAndSeed
        (
            mesh,
            electroProperties,
            initialDeltaT,
            ionicModelPtr(),
            tissuePtr->VmRef()
        );
    }

    return autoPtr<myocardiumDomainInterface>(tissuePtr.ptr());
}

} // End namespace Foam

// ************************************************************************* //
