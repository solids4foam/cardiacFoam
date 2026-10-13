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

#include "eikonalTemplateSource.H"
#include "ecgDomain.H"
#include "dimVoltage.H"
#include "mathematicalConstants.H"
#include "ionicHeterogeneity.H"
#include "IOdictionary.H"
#include "stimulusIO.H"
#include "tissueTemplates.H"

namespace Foam
{

eikonalTemplateSource::eikonalTemplateSource(const dictionary& dict)
:
    useManufacturedTemplate_(false),
    templateType_(dict.lookupOrDefault<word>("templateType", "regional")),
    weightsCalculated_(false),
    personalizedTemplatesEnabled_(false),
    personalizedTemplatesGenerated_(false)
{
    const word verifierType =
        dict.lookupOrDefault<word>("ecgVerificationModel", word::null);

    word verifierTypeModern = word::null;
    if (dict.found("verificationModel"))
    {
        verifierTypeModern = dict.subDict("verificationModel").lookupOrDefault<word>("type", "none");
    }

    if
    (
        verifierType == "manufacturedEikonalECGVerifier"
     || verifierTypeModern == "manufacturedEikonalECGVerifier"
     || verifierType == "manufacturedEikonalLeadFieldECGVerifier"
     || verifierTypeModern == "manufacturedEikonalLeadFieldECGVerifier"
    )
    {
        useManufacturedTemplate_ = true;
    }

    if (templateType_ != "regional" && templateType_ != "manufactured"
     && templateType_ != "endo" && templateType_ != "mid" && templateType_ != "epi")
    {
        FatalIOErrorInFunction(dict)
            << "templateType must be regional, manufactured, endo, mid or epi."
            << exit(FatalIOError);
    }
    if (useManufacturedTemplate_ && templateType_ != "regional"
     && templateType_ != "manufactured")
    {
        FatalIOErrorInFunction(dict)
            << "Manufactured eikonal ECG verification requires the manufactured template."
            << exit(FatalIOError);
    }
    useManufacturedTemplate_ = useManufacturedTemplate_ || templateType_ == "manufactured";
    if (useManufacturedTemplate_) templateType_ = "manufactured";

    if (dict.found("personalizedTemplates"))
    {
        if (templateType_ != "regional")
        {
            FatalIOErrorInFunction(dict)
                << "personalizedTemplates requires templateType regional."
                << exit(FatalIOError);
        }
        personalizedTemplatesDict_ = dict.subDict("personalizedTemplates");
        personalizedTemplatesEnabled_ = true;

        if (useManufacturedTemplate_)
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates cannot be combined "
                << "with a manufactured ECG configuration "
                << "(ecgVerificationModel / verificationModel.type == "
                << "manufacturedEikonalECGVerifier)."
                << exit(FatalError);
        }

        if (!personalizedTemplatesDict_.found("ionicModelConfig"))
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates requires an "
                << "ionicModelConfig sub-dictionary: an eikonal-only case "
                << "has no ionic model to inherit one from."
                << exit(FatalError);
        }

        const dictionary& ionicModelConfig =
            personalizedTemplatesDict_.subDict("ionicModelConfig");

        if (!ionicModelConfig.found("ionicModel"))
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates.ionicModelConfig "
                << "requires an 'ionicModel' entry."
                << exit(FatalError);
        }

        if (!ionicModelConfig.found("singleCellStimulus"))
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates.ionicModelConfig "
                << "requires a 'singleCellStimulus' sub-dictionary defining "
                << "the S1 pacing protocol."
                << exit(FatalError);
        }

        const label nBeats = personalizedTemplatesDict_.get<label>("nBeats");
        const scalar duration =
            personalizedTemplatesDict_.get<scalar>("duration");
        const scalar generationDt =
            personalizedTemplatesDict_.get<scalar>("dt");

        if (nBeats < 1)
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates.nBeats must be >= 1; "
                << "got " << nBeats << "."
                << exit(FatalError);
        }

        if (duration <= 0.0 || generationDt <= 0.0)
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates.duration and .dt must "
                << "both be positive; got duration=" << duration
                << ", dt=" << generationDt << "."
                << exit(FatalError);
        }

        const StimulusProtocol stim =
            stimulusIO::loadStimulusProtocol(ionicModelConfig);

        if (duration > 1e-3*stim.stimPeriodS1)
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates.duration (" << duration
                << " s) must not exceed one S1 period ("
                << 1e-3*stim.stimPeriodS1 << " s, from ionicModelConfig."
                << "singleCellStimulus.stim_period_S1=" << stim.stimPeriodS1
                << " ms): capture must stay within the final S1 beat."
                << exit(FatalError);
        }

        if (stim.nStim2 != 0)
        {
            FatalErrorInFunction
                << "eikonalECG personalizedTemplates does not support S2 "
                << "pacing (ionicModelConfig.singleCellStimulus.nstim2 must "
                << "be 0); capture must contain exactly one unambiguous "
                << "final S1 response."
                << exit(FatalError);
        }
    }
}


void eikonalTemplateSource::initialise(const ecgDomain& domain)
{
    if (templateType_ != "regional") return;
    if (personalizedTemplatesEnabled_ && !personalizedTemplatesGenerated_)
        generatePersonalizedTemplates(domain);
    if (!weightsCalculated_) calculateTransmuralWeights(domain);
}


scalar eikonalTemplateSource::manufacturedTemplateDerivative(scalar localTime) const
{
    return
        2.0*constant::mathematical::pi
       *Foam::cos(2.0*constant::mathematical::pi*localTime);
}


void eikonalTemplateSource::reconstructGradVm
(
    scalar sampleTime,
    const volScalarField& activationTime,
    const volVectorField& gradActivationTime,
    volVectorField& gradVm
) const
{
    if (&activationTime.mesh() != &gradActivationTime.mesh()
     || &activationTime.mesh() != &gradVm.mesh()
     || activationTime.dimensions() != dimTime
     || gradActivationTime.dimensions() != dimTime/dimLength
     || gradVm.dimensions() != dimVoltage/dimLength)
    {
        FatalErrorInFunction << "Template reconstruction requires matching meshes, "
            << "activation time in seconds and voltage gradient in V/m."
            << exit(FatalError);
    }
    const scalarField& activationValues = activationTime.primitiveField();
    const vectorField& gradTauValues = gradActivationTime.primitiveField();
    vectorField& gradVmValues = gradVm.primitiveFieldRef();

    forAll(gradVmValues, cellI)
    {
        const scalar localTime = sampleTime - activationValues[cellI];
        scalar dUds;

        if (useManufacturedTemplate_)
        {
            dUds = manufacturedTemplateDerivative(localTime);
        }
        else if (templateType_ != "regional")
        {
            const scalar* times = eikonalECG_templates::endoTimes;
            const scalar* values = eikonalECG_templates::endoValues;
            label count = eikonalECG_templates::numEndoSamples;
            if (templateType_ == "mid")
            {
                times = eikonalECG_templates::midTimes;
                values = eikonalECG_templates::midValues;
                count = eikonalECG_templates::numMidSamples;
            }
            else if (templateType_ == "epi")
            {
                times = eikonalECG_templates::epiTimes;
                values = eikonalECG_templates::epiValues;
                count = eikonalECG_templates::numEpiSamples;
            }
            dUds = 1e-3*eikonalECG_templates::evaluateTemplateDerivative
            (
                localTime, times, values, count
            );
        }
        else if (personalizedTemplatesEnabled_)
        {
            scalar rawDUds = 0.0;
            forAll(personalizedTemplates_, i)
            {
                rawDUds +=
                    regionWeights_[i][cellI]
                   *eikonalECG_templates::evaluateTemplateDerivative
                    (
                        localTime,
                        personalizedTemplates_[i].times.cdata(),
                        personalizedTemplates_[i].valuesMv.cdata(),
                        personalizedTemplates_[i].times.size()
                    );
            }

            // Dynamic values are mV; convert to V/s before the chain rule.
            dUds = rawDUds*1e-3;
        }
        else
        {
            // Compiled fallback is always exactly 3 fixed curves in
            // [endo, mid, epi] order, matching
            // synthesizeTransmuralBandRegions()'s own fixed ordering.
            const scalar rawDUds =
                regionWeights_[0][cellI] * eikonalECG_templates::evaluateTemplateDerivative
                (
                    localTime,
                    eikonalECG_templates::endoTimes,
                    eikonalECG_templates::endoValues,
                    eikonalECG_templates::numEndoSamples
                )
              + regionWeights_[1][cellI] * eikonalECG_templates::evaluateTemplateDerivative
                (
                    localTime,
                    eikonalECG_templates::midTimes,
                    eikonalECG_templates::midValues,
                    eikonalECG_templates::numMidSamples
                )
              + regionWeights_[2][cellI] * eikonalECG_templates::evaluateTemplateDerivative
                (
                    localTime,
                    eikonalECG_templates::epiTimes,
                    eikonalECG_templates::epiValues,
                    eikonalECG_templates::numEpiSamples
                );
            // Templates are stored in mV; field uses dimVoltage (V).
            dUds = rawDUds * 1e-3;
        }

        // Uniform-template chain rule; legacy regional weights are held fixed.
        gradVmValues[cellI] = -dUds * gradTauValues[cellI];
    }
}


const dictionary* eikonalTemplateSource::findHeterogeneityDict
(
    const dictionary& electroProperties
) const
{
    const dictionary* hetDictPtr =
        electroProperties.findDict("ionicHeterogeneity");

    if (!hetDictPtr && electroProperties.found("myocardiumSolver"))
    {
        const word myocardiumSolver(electroProperties.lookup("myocardiumSolver"));
        const dictionary* solverDictPtr =
            electroProperties.findDict(myocardiumSolver + "Coeffs");

        if (solverDictPtr)
        {
            hetDictPtr = solverDictPtr->findDict("ionicHeterogeneity");
        }
    }

    return hetDictPtr;
}


void eikonalTemplateSource::calculateTransmuralWeights(const ecgDomain& domain)
{
    const fvMesh& mesh = domain.mesh();
    weightsCalculated_ = true;

    if (useManufacturedTemplate_)
    {
        return;
    }

    IOdictionary electroProperties
    (
        IOobject
        (
            "electroProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary* hetDictPtr = findHeterogeneityDict(electroProperties);

    if (!hetDictPtr)
    {
        FatalErrorInFunction
            << "Heterogeneous templates used, but no ionicHeterogeneity "
            << "block found in electroProperties or the active myocardium "
            << "solver Coeffs dictionary."
            << exit(FatalError);
    }

    const dictionary& hetDict = *hetDictPtr;

    if (!hetDict.found("mode"))
    {
        FatalErrorInFunction
            << "eikonalECG ionicHeterogeneity has no 'mode' entry. 'mode' "
            << "is required: 'namedRegions' (always supported here) or "
            << "'cellZoneRegions' (only when personalizedTemplates is also "
            << "enabled)."
            << exit(FatalError);
    }

    const word mode(hetDict.lookup("mode"));

    const bool cellZoneRegionsSupported =
        mode == "cellZoneRegions" && personalizedTemplatesEnabled_;

    if (mode != "namedRegions" && !cellZoneRegionsSupported)
    {
        FatalErrorInFunction
            << "eikonalECG supports ionicHeterogeneity mode 'namedRegions' "
            << "always, or 'cellZoneRegions' only when personalizedTemplates "
            << "is also enabled (the fixed 3-curve compiled fallback cannot "
            << "represent an arbitrary number of regions); mode '" << mode
            << "' cannot be represented by per-region template weights here."
            << exit(FatalError);
    }

    if (mode == "cellZoneRegions")
    {
        const List<ionicHeterogeneity::NamedCellZoneRegion> czRegions =
            ionicHeterogeneity::parseNamedCellZoneRegions
            (
                hetDict.subDict("regions")
            );

        const label nCzRegions = czRegions.size();
        regionWeights_.setSize(nCzRegions);
        forAll(regionWeights_, i)
        {
            regionWeights_[i].setSize(mesh.nCells(), 0.0);
        }

        boolList claimed(mesh.nCells(), false);

        forAll(czRegions, regionIndex)
        {
            const word& zoneName = czRegions[regionIndex].cellZone;
            const label zoneId = mesh.cellZones().findZoneID(zoneName);

            if (zoneId < 0)
            {
                FatalErrorInFunction
                    << "eikonalECG personalizedTemplates: cellZone '"
                    << zoneName << "' (region '"
                    << czRegions[regionIndex].name << "') not found in "
                    << "this mesh."
                    << exit(FatalError);
            }

            const labelList& zoneCells = mesh.cellZones()[zoneId];

            forAll(zoneCells, i)
            {
                const label cellI = zoneCells[i];

                if (claimed[cellI])
                {
                    FatalErrorInFunction
                        << "eikonalECG personalizedTemplates: cell " << cellI
                        << " is claimed by more than one cellZoneRegions "
                        << "entry."
                        << exit(FatalError);
                }

                claimed[cellI] = true;
                regionWeights_[regionIndex][cellI] = 1.0;
            }
        }

        forAll(claimed, cellI)
        {
            if (!claimed[cellI])
            {
                FatalErrorInFunction
                    << "eikonalECG personalizedTemplates: cell " << cellI
                    << " is not claimed by any cellZoneRegions entry. "
                    << "Every mesh cell must belong to exactly one named "
                    << "cell zone."
                    << exit(FatalError);
            }
        }

        return;
    }

    if (!hetDict.found("transitionMode"))
    {
        FatalErrorInFunction
            << "eikonalECG ionicHeterogeneity mode namedRegions requires a "
            << "'transitionMode' entry. Supported: blend, hard."
            << exit(FatalError);
    }

    const word transitionMode(hetDict.lookup("transitionMode"));

    word smoothing;
    scalar transitionWidth = 0.0;

    if (transitionMode == "blend")
    {
        if (!hetDict.found("transitionWidth") || !hetDict.found("smoothing"))
        {
            FatalErrorInFunction
                << "eikonalECG ionicHeterogeneity mode namedRegions with "
                << "transitionMode blend requires both 'transitionWidth' "
                << "and 'smoothing' entries."
                << exit(FatalError);
        }

        transitionWidth = hetDict.get<scalar>("transitionWidth");
        smoothing = word(hetDict.lookup("smoothing"));
    }
    else if (transitionMode != "hard")
    {
        FatalErrorInFunction
            << "Unsupported eikonalECG ionicHeterogeneity transitionMode '"
            << transitionMode << "'. Supported: blend, hard."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedFieldRegion> regions =
        ionicHeterogeneity::parseNamedFieldRegions
        (
            hetDict.subDict("regions")
        );

    if (!hetDict.found("field"))
    {
        FatalErrorInFunction
            << "eikonalECG ionicHeterogeneity mode namedRegions requires a "
            << "'field' entry."
            << exit(FatalError);
    }

    const word fieldName(hetDict.lookup("field"));
    autoPtr<volScalarField> tReadPtr;
    const volScalarField* tPtr = mesh.cfindObject<volScalarField>(fieldName);

    if (!tPtr)
    {
        // Transmural-distance fields are static anatomical inputs stored in
        // the latest time directory that contains the field.
        const word fieldInstance =
            mesh.time().findInstance(fileName::null, fieldName);

        tReadPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    fieldInstance,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        tPtr = &tReadPtr();
    }

    if (!tPtr)
    {
        FatalErrorInFunction
            << "Could not find transmural distance field '" << fieldName << "'."
            << exit(FatalError);
    }

    const scalarField& tField = tPtr->primitiveField();

    const label nRegions = regions.size();
    regionWeights_.setSize(nRegions);
    forAll(regionWeights_, i)
    {
        regionWeights_[i].setSize(mesh.nCells(), 0.0);
    }

    // Small, fixed-size lookup built once, not per-cell.
    List<word> regionNames(nRegions);
    forAll(regions, i)
    {
        regionNames[i] = regions[i].name;
    }

    forAll(tField, cellI)
    {
        const List<ionicHeterogeneity::NamedRegionWeight> sparse =
            ionicHeterogeneity::namedRegionWeightsAt
            (
                tField[cellI], regions, transitionWidth, smoothing, transitionMode
            );

        forAll(sparse, s)
        {
            forAll(regionNames, i)
            {
                if (regionNames[i] == sparse[s].name)
                {
                    regionWeights_[i][cellI] = sparse[s].weight;
                    break;
                }
            }
        }
    }
}


void eikonalTemplateSource::generatePersonalizedTemplates(const ecgDomain& domain)
{
    const fvMesh& mesh = domain.mesh();

    IOdictionary electroProperties
    (
        IOobject
        (
            "electroProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary* hetDictPtr = findHeterogeneityDict(electroProperties);

    if (!hetDictPtr)
    {
        FatalErrorInFunction
            << "eikonalECG personalizedTemplates requires an "
            << "ionicHeterogeneity block in electroProperties or the "
            << "active myocardium solver Coeffs dictionary, but none was "
            << "found."
            << exit(FatalError);
    }

    const dictionary& hetDict = *hetDictPtr;

    if (!hetDict.found("mode"))
    {
        FatalErrorInFunction
            << "eikonalECG personalizedTemplates requires ionicHeterogeneity "
            << "'mode' to be set: 'namedRegions' or 'cellZoneRegions'."
            << exit(FatalError);
    }

    const word mode(hetDict.lookup("mode"));

    if (mode != "namedRegions" && mode != "cellZoneRegions")
    {
        FatalErrorInFunction
            << "eikonalECG personalizedTemplates supports ionicHeterogeneity "
            << "mode 'namedRegions' or 'cellZoneRegions' only; mode '"
            << mode << "' cannot be used to generate anchor templates."
            << exit(FatalError);
    }

    const dictionary& ionicModelConfig =
        personalizedTemplatesDict_.subDict("ionicModelConfig");

    const label nBeats = personalizedTemplatesDict_.get<label>("nBeats");
    const scalar duration = personalizedTemplatesDict_.get<scalar>("duration");
    const scalar generationDt = personalizedTemplatesDict_.get<scalar>("dt");

    personalizedTemplates_ =
        eikonalECG_templates::generatePersonalizedTemplates
        (
            ionicModelConfig,
            hetDict,
            nBeats,
            duration,
            generationDt
        );

    personalizedTemplatesGenerated_ = true;
}


} // End namespace Foam
