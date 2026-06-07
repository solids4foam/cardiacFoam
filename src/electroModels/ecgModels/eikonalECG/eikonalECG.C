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
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicHeterogeneity.H"
#include "IOdictionary.H"
#include "tissueTemplates.H"

namespace Foam
{

defineTypeNameWithName(eikonalECG, "eikonalECG");
defineDebugSwitch(eikonalECG, 0);
addToRunTimeSelectionTable(ecgSolver, eikonalECG, dictionary);


eikonalECG::eikonalECG(const dictionary& dict)
:
    useManufacturedTemplate_(false),
    startTime_(0.0),
    endTime_(0.0),
    deltaT_(0.0),
    report_(dict.lookupOrDefault<Switch>("report", true)),
    written_(false),
    lastValues_(),
    outputPtr_(),
    VmPtr_(),
    leadVectorsCalculated_(false),
    leadVectors_(),
    weightsCalculated_(false)
{
    const word verifierType =
        dict.lookupOrDefault<word>("ecgVerificationModel", word::null);

    if
    (
        dict.found("manufacturedEikonalECG")
     || verifierType == "eikonalECGManufacturedVerifier"
    )
    {
        useManufacturedTemplate_ = true;
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
}


scalar eikonalECG::manufacturedTemplateValue(scalar localTime) const
{
    return Foam::sin(2.0*constant::mathematical::pi*localTime);
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
            VmValues[cellI] = manufacturedTemplateValue(localTime);
        }
        else
        {
            if (wEndo_[cellI] > 0.5)
            {
                VmValues[cellI] = eikonalECG_templates::evaluateTemplate
                (
                    localTime,
                    eikonalECG_templates::endoTimes,
                    eikonalECG_templates::endoValues,
                    eikonalECG_templates::numEndoSamples
                );
            }
            else if (wMid_[cellI] > 0.5)
            {
                VmValues[cellI] = eikonalECG_templates::evaluateTemplate
                (
                    localTime,
                    eikonalECG_templates::midTimes,
                    eikonalECG_templates::midValues,
                    eikonalECG_templates::numMidSamples
                );
            }
            else
            {
                VmValues[cellI] = eikonalECG_templates::evaluateTemplate
                (
                    localTime,
                    eikonalECG_templates::epiTimes,
                    eikonalECG_templates::epiValues,
                    eikonalECG_templates::numEpiSamples
                );
            }
        }
    }

    Vm.correctBoundaryConditions();
}


void eikonalECG::calculatePseudoECG
(
    const ecgDomain& domain,
    const volScalarField& Vm,
    scalarField& values
) const
{
    const label nElectrodes = domain.electrodePositions().size();

    const tmp<volVectorField> tgradVm = fvc::grad(Vm);
    const vectorField& gradVm = tgradVm().primitiveField();

    values.setSize(nElectrodes);
    values = 0.0;

    for (label electrodeI = 0; electrodeI < nElectrodes; ++electrodeI)
    {
        const List<vector>& z = leadVectors_[electrodeI];
        scalar ecgVal = 0.0;

        forAll(gradVm, cellI)
        {
            ecgVal += gradVm[cellI] & z[cellI];
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

    wEndo_.setSize(mesh.nCells(), 1.0);
    wMid_.setSize(mesh.nCells(), 0.0);
    wEpi_.setSize(mesh.nCells(), 0.0);
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

    if (!electroProperties.found("ionicHeterogeneity"))
    {
        FatalErrorInFunction
            << "Heterogeneous templates used, but no ionicHeterogeneity "
            << "block found in electroProperties."
            << exit(FatalError);
    }

    const dictionary& hetDict = electroProperties.subDict("ionicHeterogeneity");
    const word mode = hetDict.lookupOrDefault<word>("mode", "transmuralBands");
    const word transitionMode = hetDict.lookupOrDefault<word>("transitionMode", "blend");

    if (transitionMode != "hard")
    {
        FatalErrorInFunction
            << "eikonalECG transmural heterogeneity requires transitionMode 'hard' "
            << "(found '" << transitionMode << "') because we are stepping between 3 distinct voltage curves."
            << exit(FatalError);
    }

    const word fieldName = hetDict.lookupOrDefault<word>("field", "t");
    const volScalarField* tPtr = mesh.cfindObject<volScalarField>(fieldName);

    if (!tPtr)
    {
        FatalErrorInFunction
            << "Could not find transmural distance field '" << fieldName << "'."
            << exit(FatalError);
    }

    const scalar endoMInterface = hetDict.lookupOrDefault<scalar>("endoMInterface", 0.3);
    const scalar mEpiInterface = hetDict.lookupOrDefault<scalar>("mEpiInterface", 0.7);
    const scalar transitionWidth = hetDict.lookupOrDefault<scalar>("transitionWidth", 0.1);
    const word smoothing = hetDict.lookupOrDefault<word>("smoothing", "smoothstep");

    ionicHeterogeneity::validateTransmuralBandConfig
    (
        endoMInterface, mEpiInterface, transitionWidth, smoothing, transitionMode
    );

    const scalarField& tField = tPtr->primitiveField();

    forAll(tField, cellI)
    {
        const ionicHeterogeneity::TransmuralBandWeights w =
            ionicHeterogeneity::transmuralBandWeights
            (
                tField[cellI], endoMInterface, mEpiInterface,
                transitionWidth, smoothing, transitionMode
            );

        wEndo_[cellI] = w.endo;
        wMid_[cellI]  = w.mCell;
        wEpi_[cellI]  = w.epi;
    }
}


volScalarField& eikonalECG::surrogateVm(ecgDomain& domain)
{
    if (!VmPtr_.valid())
    {
        VmPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "eikonalECG_Vm",
                    domain.mesh().time().timeName(),
                    domain.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                domain.mesh(),
                dimensionedScalar
                (
                    "zero",
                    dimVoltage,
                    0.0
                ),
                "zeroGradient"
            )
        );
    }

    return VmPtr_.ref();
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

    if (written_)
    {
        values = lastValues_;
        return;
    }

    if (!leadVectorsCalculated_)
    {
        calculateLeadVectors(domain);
    }

    if (!weightsCalculated_)
    {
        calculateTransmuralWeights(domain);
    }

    if (!outputPtr_.valid())
    {
        outputPtr_ =
            ecgModelIO::openTimeSeries
            (
                domain.mesh().time().globalPath()/"postProcessing",
                "eikonalECG.dat",
                domain.electrodeNames()
            );
    }

    const volScalarField& activationTime = domain.activationTime();
    volScalarField& Vm = surrogateVm(domain);

    scalar sampleTime = startTime_;
    label sampleI = 0;

    while (sampleTime <= endTime_ + SMALL)
    {
        reconstructVm(sampleTime, activationTime, Vm);
        calculatePseudoECG(domain, Vm, lastValues_);
        ecgModelIO::writeRow(outputPtr_.ref(), sampleTime, lastValues_);
        domain.recordVerification(sampleTime, lastValues_);

        ++sampleI;
        sampleTime = startTime_ + sampleI*deltaT_;
    }

    values = lastValues_;
    written_ = true;

    if (report_)
    {
        Info<< "eikonalECG: wrote sampled ECG from t=" << startTime_
            << " to t=" << endTime_ << " with deltaT=" << deltaT_
            << " to postProcessing/eikonalECG.dat" << nl << endl;

        if (useManufacturedTemplate_)
        {
            Info<< "eikonalECG: using fixed manufactured eikonal ECG "
                << "template because no template dictionary was supplied."
                << nl << endl;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
