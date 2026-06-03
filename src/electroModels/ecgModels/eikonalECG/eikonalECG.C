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
#include "eikonalVerification/manufacturedEikonalReference.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameWithName(eikonalECG, "eikonalECG");
defineDebugSwitch(eikonalECG, 0);
addToRunTimeSelectionTable(ecgSolver, eikonalECG, dictionary);


eikonalECG::eikonalECG(const dictionary& dict)
:
    templateTimes_(),
    templateValues_(),
    useManufacturedTemplate_(false),
    startTime_(0.0),
    endTime_(0.0),
    deltaT_(0.0),
    report_(dict.lookupOrDefault<Switch>("report", true)),
    written_(false),
    lastValues_(),
    outputPtr_(),
    VmPtr_()
{
    if (dict.found("template"))
    {
        const dictionary& templateDict = dict.subDict("template");
        templateTimes_ = templateDict.get<scalarField>("times");
        templateValues_ = templateDict.get<scalarField>("values");
        validateTemplate();
    }
    else
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
        else
        {
            FatalErrorInFunction
                << "eikonalECG requires a template dictionary unless "
                << "manufacturedEikonalECG verification is enabled."
                << exit(FatalError);
        }
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


void eikonalECG::validateTemplate() const
{
    if (useManufacturedTemplate_)
    {
        return;
    }

    if (templateTimes_.size() != templateValues_.size())
    {
        FatalErrorInFunction
            << "eikonalECG template.times and template.values must have "
            << "the same size."
            << exit(FatalError);
    }

    if (templateTimes_.size() < 2)
    {
        FatalErrorInFunction
            << "eikonalECG template requires at least two samples."
            << exit(FatalError);
    }

    for (label i = 1; i < templateTimes_.size(); ++i)
    {
        if (templateTimes_[i] <= templateTimes_[i - 1])
        {
            FatalErrorInFunction
                << "eikonalECG template.times must be strictly increasing."
                << exit(FatalError);
        }
    }
}


scalar eikonalECG::templateValue(scalar localTime) const
{
    if (useManufacturedTemplate_)
    {
        return manufacturedEikonalTemplateValue(localTime);
    }

    if (localTime <= templateTimes_.first())
    {
        return templateValues_.first();
    }

    if (localTime >= templateTimes_.last())
    {
        return templateValues_.last();
    }

    label lo = 0;
    label hi = templateTimes_.size() - 1;

    while (hi - lo > 1)
    {
        const label mid = (lo + hi)/2;

        if (templateTimes_[mid] <= localTime)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    const scalar t0 = templateTimes_[lo];
    const scalar t1 = templateTimes_[hi];
    const scalar alpha = (localTime - t0)/(t1 - t0);

    return (1.0 - alpha)*templateValues_[lo] + alpha*templateValues_[hi];
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
        VmValues[cellI] =
            templateValue(sampleTime - activationValues[cellI]);
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
    const fvMesh& mesh = domain.mesh();
    const List<vector>& electrodePositions = domain.electrodePositions();

    const tmp<volVectorField> tgradVm = fvc::grad(Vm);
    const vectorField& gradVm = tgradVm().primitiveField();

    const scalarField& volumes = mesh.V();
    const vectorField& cellCentres = mesh.C().primitiveField();
    const tensorField& conductivityField =
        domain.conductivity().primitiveField();

    const label nElectrodes = electrodePositions.size();

    values.setSize(nElectrodes);
    values = 0.0;

    forAll(cellCentres, cellI)
    {
        const vector dipole =
            (conductivityField[cellI] & gradVm[cellI]) * volumes[cellI];

        for (label electrodeI = 0; electrodeI < nElectrodes; ++electrodeI)
        {
            const vector rVec =
                cellCentres[cellI] - electrodePositions[electrodeI];
            const scalar r = mag(rVec);

            if (r > VSMALL)
            {
                values[electrodeI] += (dipole & rVec)/(r*r*r);
            }
        }
    }

    for (label electrodeI = 0; electrodeI < nElectrodes; ++electrodeI)
    {
        reduce(values[electrodeI], sumOp<scalar>());
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
