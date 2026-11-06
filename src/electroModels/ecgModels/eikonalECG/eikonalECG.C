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

#include "eikonalECG.H"

#include "dimVoltage.H"
#include "ecgDomain.H"
#include "ecgModelIO.H"
#include "mathematicalConstants.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "ionicHeterogeneity.H"
#include "tissueTemplates.H"
#include "manufacturedTemplate.H"

namespace Foam
{

defineTypeNameWithName(eikonalECG, "eikonalECG");
defineDebugSwitch(eikonalECG, 0);
addToRunTimeSelectionTable(ecgSolver, eikonalECG, dictionary);


eikonalECG::eikonalECG(const dictionary& dict)
:
    useManufacturedTemplate_(false),
    manufacturedSineTemplate_(false),
    manufacturedTemplateDuration_(0),
    manufacturedTemplateRise_(0),
    startTime_(0.0),
    endTime_(0.0),
    deltaT_(0.0),
    sigmaE_(dict.lookupOrDefault<scalar>("sigmaExtracellular", 0.0)),
    lastValues_(),
    outputPtr_(),
    leadVectorsCalculated_(false),
    leadVectors_(),
    weightsCalculated_(false),
    regionWeights_(),
    regionWeightGradients_(),
    compiledCurveIndex_(),
    checkGradVm_(dict.lookupOrDefault<Switch>("checkGradVm", false)),
    personalizedTemplatesEnabled_(false),
    personalizedTemplatesGenerated_(false)
{
    const word verifierType =
        dict.lookupOrDefault<word>("ecgVerificationModel", word::null);

    word verifierTypeModern = word::null;
    if (dict.found("verificationModel"))
    {
        verifierTypeModern =
            dict.subDict("verificationModel").lookupOrDefault<word>
            (
                "type",
                "none"
            );
    }

    if
    (
        verifierType == "manufacturedEikonalECGVerifier"
     || verifierTypeModern == "manufacturedEikonalECGVerifier"
    )
    {
        useManufacturedTemplate_ = true;
        eikonalECG_templates::readManufacturedTemplate
        (
            dict.subOrEmptyDict("verificationModel"),
            manufacturedSineTemplate_,
            manufacturedTemplateDuration_,
            manufacturedTemplateRise_
        );
    }

    const dictionary& samplingDict = dict.subDict("sampling");
    startTime_ = samplingDict.get<scalar>("start");
    endTime_ = samplingDict.get<scalar>("end");
    deltaT_ = samplingDict.get<scalar>("deltaT");

    if (deltaT_ <= 0.0)
    {
        FatalErrorInFunction
            << "eikonalECG sampling.deltaT must be positive."
            << exit(FatalError);
    }

    if (endTime_ < startTime_)
    {
        FatalErrorInFunction
            << "eikonalECG sampling.end must be greater than or equal to "
            << "sampling.start."
            << exit(FatalError);
    }

    if (dict.found("personalizedTemplates"))
    {
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

        ionicTemplates::checkControls(personalizedTemplatesDict_);
    }
}


void eikonalECG::evaluateRegionTemplate
(
    label regionI,
    scalar localTime,
    scalar& value,
    scalar& slope
) const
{
    // Templates are stored in mV; fields use dimVoltage (V).
    if (personalizedTemplatesEnabled_)
    {
        const ionicTemplates::DynamicTemplate& tpl =
            personalizedTemplates_[regionI];
        const label n = tpl.times.size();
        value = 1e-3*ionicTemplates::evaluateTemplate
        (
            localTime, tpl.times.cdata(), tpl.valuesMv.cdata(), n
        );
        slope = 1e-3*eikonalECG_templates::evaluateTemplateDerivative
        (
            localTime, tpl.times.cdata(), tpl.valuesMv.cdata(), n
        );
        return;
    }

    // Compiled fallback: 3 fixed curves [endo, mid, epi], selected by
    // the region's name (compiledCurveIndex_), not by its position.
    const scalar* times = nullptr;
    const scalar* values = nullptr;
    int n = 0;

    switch (compiledCurveIndex_[regionI])
    {
        case 0:
            times = eikonalECG_templates::endoTimes;
            values = eikonalECG_templates::endoValues;
            n = eikonalECG_templates::numEndoSamples;
            break;
        case 1:
            times = eikonalECG_templates::midTimes;
            values = eikonalECG_templates::midValues;
            n = eikonalECG_templates::numMidSamples;
            break;
        case 2:
            times = eikonalECG_templates::epiTimes;
            values = eikonalECG_templates::epiValues;
            n = eikonalECG_templates::numEpiSamples;
            break;
        default:
            FatalErrorInFunction
                << "Compiled eikonalECG templates have 3 regions "
                << "[endo, mid, epi]; region " << regionI << " requested."
                << exit(FatalError);
    }

    value = 1e-3*ionicTemplates::evaluateTemplate
    (
        localTime, times, values, n
    );
    slope = 1e-3*eikonalECG_templates::evaluateTemplateDerivative
    (
        localTime, times, values, n
    );
}


void eikonalECG::calculateRegionWeightGradients(const fvMesh& mesh)
{
    regionWeightGradients_.setSize(regionWeights_.size());

    volScalarField w
    (
        IOobject
        (
            "eikonalECG_regionWeight",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(regionWeights_, i)
    {
        w.primitiveFieldRef() = regionWeights_[i];
        w.correctBoundaryConditions();
        regionWeightGradients_[i] = fvc::grad(w)().primitiveField();
    }
}


void eikonalECG::reconstructGradVm
(
    scalar sampleTime,
    const volScalarField& activationTime,
    const volVectorField& gradActivationTime,
    volVectorField& gradVm
) const
{
    const scalarField& activationValues = activationTime.primitiveField();
    const vectorField& gradTauValues = gradActivationTime.primitiveField();
    vectorField& gradVmValues = gradVm.primitiveFieldRef();

    forAll(gradVmValues, cellI)
    {
        const scalar localTime = sampleTime - activationValues[cellI];

        if (useManufacturedTemplate_)
        {
            const scalar dUds =
                eikonalECG_templates::evaluateManufacturedTemplateDerivative
                (
                    localTime,
                    manufacturedSineTemplate_,
                    manufacturedTemplateDuration_,
                    manufacturedTemplateRise_
                );

            // Vm(x,t) = U(t - tau(x))  =>  gradVm = -dU/ds(t-tau) * gradTau
            gradVmValues[cellI] = -dUds*gradTauValues[cellI];
            continue;
        }

        // Vm(x,t) = sum_k w_k(x) U_k(t - tau(x))
        //   => gradVm = -sum_k w_k dU_k/ds gradTau + sum_k U_k grad(w_k)
        scalar weightedSlope = 0.0;
        vector weightTerm = Zero;

        forAll(regionWeights_, i)
        {
            scalar value = 0.0;
            scalar slope = 0.0;
            evaluateRegionTemplate(i, localTime, value, slope);

            weightedSlope += regionWeights_[i][cellI]*slope;
            weightTerm += value*regionWeightGradients_[i][cellI];
        }

        gradVmValues[cellI] = -weightedSlope*gradTauValues[cellI] + weightTerm;
    }
}


void eikonalECG::reconstructVm
(
    scalar sampleTime,
    const volScalarField& activationTime,
    volScalarField& Vm
) const
{
    const scalarField& activationValues = activationTime.primitiveField();
    scalarField& VmValues = Vm.primitiveFieldRef();

    forAll(VmValues, cellI)
    {
        const scalar localTime = sampleTime - activationValues[cellI];

        if (useManufacturedTemplate_)
        {
            VmValues[cellI] = eikonalECG_templates::evaluateManufacturedTemplate
            (
                localTime,
                manufacturedSineTemplate_,
                manufacturedTemplateDuration_,
                manufacturedTemplateRise_
            );
            continue;
        }

        scalar v = 0.0;
        forAll(regionWeights_, i)
        {
            scalar value = 0.0;
            scalar slope = 0.0;
            evaluateRegionTemplate(i, localTime, value, slope);
            v += regionWeights_[i][cellI]*value;
        }
        VmValues[cellI] = v;
    }

    Vm.correctBoundaryConditions();
}


void eikonalECG::calculatePseudoECG
(
    const ecgDomain& domain,
    const volVectorField& gradVm,
    scalarField& values
) const
{
    const label nElectrodes = domain.electrodePositions().size();
    const vectorField& gradVmValues = gradVm.primitiveField();

    values.setSize(nElectrodes);
    values = 0.0;

    for (label electrodeI = 0; electrodeI < nElectrodes; ++electrodeI)
    {
        const List<vector>& z = leadVectors_[electrodeI];
        scalar ecgVal = 0.0;

        forAll(gradVmValues, cellI)
        {
            ecgVal += gradVmValues[cellI] & z[cellI];
        }

        reduce(ecgVal, sumOp<scalar>());

        values[electrodeI] = ecgVal;
    }
}

void eikonalECG::calculateLeadVectors(const ecgDomain& domain)
{
    const fvMesh& mesh = domain.mesh();
    const List<vector>& electrodePositions = domain.electrodePositions();
    const label nElectrodes = electrodePositions.size();

    const scalarField& volumes = mesh.V();
    const vectorField& cellCentres = mesh.C().primitiveField();
    const tensorField& conductivityField =
        domain.conductivity().primitiveField();

    leadVectors_.setSize(nElectrodes);
    forAll(leadVectors_, i)
    {
        leadVectors_[i].setSize(cellCentres.size(), vector::zero);
    }

    forAll(cellCentres, cellI)
    {
        for (label electrodeI = 0; electrodeI < nElectrodes; ++electrodeI)
        {
            const vector rVec =
                cellCentres[cellI] - electrodePositions[electrodeI];
            const scalar r = mag(rVec);

            if (r > VSMALL)
            {
                leadVectors_[electrodeI][cellI] =
                    (conductivityField[cellI] & rVec)*(volumes[cellI]/(r*r*r));
            }
        }
    }

    leadVectorsCalculated_ = true;
}

void eikonalECG::calculateTransmuralWeights(const ecgDomain& domain)
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

    const dictionary* hetDictPtr =
        ionicTemplates::findHeterogeneityDict(electroProperties);

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

    const word mode(hetDict.get<word>("mode"));

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

    regionWeights_ = ionicTemplates::regionWeights(mesh, hetDict);

    // Compiled-template path: map each region to one of the 3 compiled
    // curves [endo, mid, epi] by name, in the order regionWeights uses
    // (ionicHeterogeneity::parseNamedFieldRegions), so a subset such as
    // endocardial and epicardial only selects the right curves.
    if (!personalizedTemplatesEnabled_)
    {
        const List<ionicHeterogeneity::NamedFieldRegion> regions =
            ionicHeterogeneity::parseNamedFieldRegions
            (
                hetDict.subDict("regions")
            );

        compiledCurveIndex_.setSize(regions.size());
        forAll(regions, i)
        {
            const word& name = regions[i].name;
            if (name == "endocardialCells")
            {
                compiledCurveIndex_[i] = 0;
            }
            else if (name == "mCells")
            {
                compiledCurveIndex_[i] = 1;
            }
            else if (name == "epicardialCells")
            {
                compiledCurveIndex_[i] = 2;
            }
            else
            {
                FatalErrorInFunction
                    << "eikonalECG compiled templates exist for regions "
                    << "endocardialCells, mCells and epicardialCells; region '"
                    << name << "' has no compiled template. Use "
                    << "personalizedTemplates for other regions."
                    << exit(FatalError);
            }
        }
    }
}


void eikonalECG::generatePersonalizedTemplates(const ecgDomain& domain)
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

    const dictionary* hetDictPtr =
        ionicTemplates::findHeterogeneityDict(electroProperties);

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

    const word mode(hetDict.get<word>("mode"));

    if (mode != "namedRegions" && mode != "cellZoneRegions")
    {
        FatalErrorInFunction
            << "eikonalECG personalizedTemplates supports ionicHeterogeneity "
            << "mode 'namedRegions' or 'cellZoneRegions' only; mode '"
            << mode << "' cannot be used to generate anchor templates."
            << exit(FatalError);
    }

    personalizedTemplates_ =
        ionicTemplates::generate
        (
            personalizedTemplatesDict_,
            hetDict,
            mesh.time()
        );

    personalizedTemplatesGenerated_ = true;
}


void eikonalECG::solve
(
    ecgDomain& domain,
    scalar t0,
    scalar dt,
    scalarField& values
)
{
    (void)t0;
    (void)dt;

    if (outputPtr_.valid())
    {
        values = lastValues_;
        return;
    }

    if (!leadVectorsCalculated_)
    {
        calculateLeadVectors(domain);
    }

    if (personalizedTemplatesEnabled_ && !personalizedTemplatesGenerated_)
    {
        generatePersonalizedTemplates(domain);
    }

    if (!weightsCalculated_)
    {
        calculateTransmuralWeights(domain);
    }

    if
    (
        !useManufacturedTemplate_
     && regionWeightGradients_.size() != regionWeights_.size()
    )
    {
        calculateRegionWeightGradients(domain.mesh());
    }

    outputPtr_ =
        ecgModelIO::openTimeSeries
        (
            domain.mesh().time().globalPath()/"postProcessing",
            "eikonalECG.dat",
            domain.electrodeNames()
        );

    const volScalarField& activationTime = domain.activationTime();

    // activationTime is fixed for the whole sampling loop below (the eikonal
    // solve is steady-state; only sampleTime varies), so its gradient only
    // needs computing once. Reusing it here (rather than reconstructing Vm
    // pointwise and differencing it afresh each sample via fvc::grad(Vm)) is
    // also what recovers the full accuracy of the underlying activation-time
    // solve in the ECG output -- see reconstructGradVm.
    const volVectorField gradActivationTime(fvc::grad(activationTime));
    volVectorField gradVm
    (
        IOobject
        (
            "eikonalECG_gradVm",
            domain.mesh().time().timeName(),
            domain.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        domain.mesh(),
        dimensionedVector("zero", dimVoltage/dimLength, Zero)
    );

    // Diagnostic only: Vm reconstructed on the mesh and differentiated by
    // fvc::grad, as an independent check of the product-rule gradient.
    autoPtr<volScalarField> checkVmPtr;
    scalarField checkValues;
    scalar checkMaxDiff = 0.0;
    scalar checkMaxSignal = 0.0;
    autoPtr<OFstream> checkOutputPtr;

    if (checkGradVm_)
    {
        checkOutputPtr =
            ecgModelIO::openTimeSeries
            (
                domain.mesh().time().globalPath()/"postProcessing",
                "eikonalECG_checkGradVm.dat",
                domain.electrodeNames()
            );

        checkVmPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "eikonalECG_checkVm",
                    domain.mesh().time().timeName(),
                    domain.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                domain.mesh(),
                dimensionedScalar("zero", dimVoltage, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    scalar sampleTime = startTime_;
    label sampleI = 0;

    while (sampleTime <= endTime_ + SMALL)
    {
        reconstructGradVm(sampleTime, activationTime, gradActivationTime, gradVm);
        calculatePseudoECG(domain, gradVm, lastValues_);

        if (checkGradVm_)
        {
            reconstructVm(sampleTime, activationTime, checkVmPtr.ref());
            const volVectorField checkGradVm(fvc::grad(checkVmPtr()));
            calculatePseudoECG(domain, checkGradVm, checkValues);

            scalarField scaledCheck(checkValues);
            if (sigmaE_ > VSMALL)
            {
                scaledCheck *= 1.0/(4.0*constant::mathematical::pi*sigmaE_);
            }
            ecgModelIO::writeRow(checkOutputPtr.ref(), sampleTime, scaledCheck);

            forAll(lastValues_, eI)
            {
                checkMaxDiff =
                    max(checkMaxDiff, mag(lastValues_[eI] - checkValues[eI]));
                checkMaxSignal = max(checkMaxSignal, mag(lastValues_[eI]));
            }
        }

        if (sigmaE_ > VSMALL)
        {
            const scalar norm =
                1.0 / (4.0 * constant::mathematical::pi * sigmaE_);
            forAll(lastValues_, eI)
            {
                lastValues_[eI] *= norm;
            }
        }

        ecgModelIO::writeRow(outputPtr_.ref(), sampleTime, lastValues_);
        domain.recordVerification(sampleTime, lastValues_);

        ++sampleI;
        sampleTime = startTime_ + sampleI*deltaT_;
    }

    values = lastValues_;

    if (checkGradVm_)
    {
        // Same normalisation for both, so the ratio is unit-free.
        Info<< "eikonalECG checkGradVm: max |ECG(product rule) - "
            << "ECG(fvc::grad Vm)| / max |ECG| = "
            << checkMaxDiff/max(checkMaxSignal, VSMALL)
            << " (before 1/(4 pi sigma) scaling)" << nl << endl;
    }

    Info<< "eikonalECG: wrote sampled ECG from t=" << startTime_
        << " to t=" << endTime_ << " with deltaT=" << deltaT_
        << " to postProcessing/eikonalECG.dat" << nl << endl;

    if (useManufacturedTemplate_)
    {
        Info<< "eikonalECG: using the manufactured eikonal ECG template "
            << "selected by the ECG verificationModel."
            << nl << endl;
    }
}

} // End namespace Foam

// ************************************************************************* //
