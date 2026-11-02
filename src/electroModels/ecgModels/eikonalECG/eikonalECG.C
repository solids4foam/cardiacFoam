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
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
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
    lastValues_(),
    outputPtr_(),
    leadField_(dict),
    weightsCalculated_(false),
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
        scalar dUds;

        if (useManufacturedTemplate_)
        {
            dUds = eikonalECG_templates::evaluateManufacturedTemplateDerivative
            (
                localTime,
                manufacturedSineTemplate_,
                manufacturedTemplateDuration_,
                manufacturedTemplateRise_
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

        // Vm(x,t) = U(t - tau(x))  =>  gradVm = -dU/ds(t-tau) * gradTau
        gradVmValues[cellI] = -dUds * gradTauValues[cellI];
    }
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

    if (personalizedTemplatesEnabled_ && !personalizedTemplatesGenerated_)
    {
        generatePersonalizedTemplates(domain);
    }

    if (!weightsCalculated_)
    {
        calculateTransmuralWeights(domain);
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

    scalar sampleTime = startTime_;
    label sampleI = 0;

    while (sampleTime <= endTime_ + SMALL)
    {
        reconstructGradVm(sampleTime, activationTime, gradActivationTime, gradVm);
        leadField_.evaluate(domain, gradVm.primitiveField(), lastValues_);

        ecgModelIO::writeRow(outputPtr_.ref(), sampleTime, lastValues_);
        domain.recordVerification(sampleTime, lastValues_);

        ++sampleI;
        sampleTime = startTime_ + sampleI*deltaT_;
    }

    values = lastValues_;

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
