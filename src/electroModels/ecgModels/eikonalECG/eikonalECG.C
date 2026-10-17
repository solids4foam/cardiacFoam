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
#include <cmath>

namespace Foam
{
defineTypeNameWithName(eikonalECG, "eikonalECG");
defineDebugSwitch(eikonalECG, 0);
addToRunTimeSelectionTable(ecgSolver, eikonalECG, dictionary);

eikonalECG::eikonalECG(const dictionary& dict)
:
    templates_(dict),
    startTime_(dict.subDict("sampling").get<scalar>("start")),
    endTime_(dict.subDict("sampling").get<scalar>("end")),
    deltaT_(dict.subDict("sampling").get<scalar>("deltaT")),
    sigmaE_(dict.lookupOrDefault<scalar>("sigmaExtracellular", 0.0)),
    lastValues_(), outputPtr_(), sampled_(false),
    leadVectorsCalculated_(false), leadVectors_()
{
    if (!std::isfinite(startTime_) || !std::isfinite(endTime_)
     || !std::isfinite(deltaT_) || deltaT_ <= 0 || endTime_ < startTime_)
    {
        FatalIOErrorInFunction(dict) << "Invalid eikonalECG sampling interval."
            << exit(FatalIOError);
    }
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

    if (sampled_)
    {
        values = lastValues_;
        return;
    }

    if (!leadVectorsCalculated_)
    {
        calculateLeadVectors(domain);
    }

    templates_.initialise(domain);

    if (Pstream::master()) outputPtr_ =
        ecgModelIO::openTimeSeries
        (
            domain.mesh().time().globalPath()/"postProcessing",
            "eikonalECG.dat",
            domain.electrodeNames()
        );

    const volScalarField& activationTime = domain.activationTime();

    // The steady activation field is fixed throughout the ECG sampling loop.
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
        templates_.reconstructGradVm(sampleTime, activationTime, gradActivationTime, gradVm);
        calculatePseudoECG(domain, gradVm, lastValues_);

        if (sigmaE_ > VSMALL)
        {
            const scalar norm =
                1.0 / (4.0 * constant::mathematical::pi * sigmaE_);
            forAll(lastValues_, eI)
            {
                lastValues_[eI] *= norm;
            }
        }

        if (Pstream::master()) ecgModelIO::writeRow(outputPtr_.ref(), sampleTime, lastValues_);
        domain.recordVerification(sampleTime, lastValues_);

        ++sampleI;
        sampleTime = startTime_ + sampleI*deltaT_;
    }

    sampled_ = true;
    values = lastValues_;

    Info<< "eikonalECG: wrote sampled ECG from t=" << startTime_
        << " to t=" << endTime_ << " with deltaT=" << deltaT_
        << " to postProcessing/eikonalECG.dat" << nl << endl;

    if (templates_.templateType() == "manufactured")
    {
        Info<< "eikonalECG: using fixed manufactured eikonal ECG "
            << "template."
            << nl << endl;
    }
}

} // End namespace Foam

// ************************************************************************* //
