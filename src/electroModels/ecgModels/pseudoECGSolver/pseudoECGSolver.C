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

#include "pseudoECGSolver.H"

#include "ecgDomain.H"
#include "ecgModelIO.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameWithName(pseudoECGSolver, "pseudoECG");
defineDebugSwitch(pseudoECGSolver, 0);
addToRunTimeSelectionTable(ecgSolver, pseudoECGSolver, dictionary);


pseudoECGSolver::pseudoECGSolver(const dictionary& dict)
:
    sigmaE_(dict.lookupOrDefault<scalar>("sigmaExtracellular", 0.0)),
    reportedConductivitySource_(false),
    leadVectorsCalculated_(false),
    hasOwnSampling_(false),
    startTime_(0.0),
    endTime_(0.0),
    deltaT_(0.0),
    nextSampleTime_(0.0),
    outputPtr_()
{
    if (const dictionary* samplingDictPtr = dict.findDict("sampling"))
    {
        hasOwnSampling_ = true;

        const dictionary& samplingDict = *samplingDictPtr;
        startTime_ = samplingDict.get<scalar>("start");
        endTime_ = samplingDict.get<scalar>("end");
        deltaT_ = samplingDict.get<scalar>("deltaT");

        if (deltaT_ <= 0.0)
        {
            FatalErrorInFunction
                << "pseudoECG sampling.deltaT must be positive."
                << exit(FatalError);
        }

        nextSampleTime_ = startTime_;
    }
}


void pseudoECGSolver::calculateLeadVectors(const ecgDomain& domain)
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


void pseudoECGSolver::solve
(
    ecgDomain& domain,
    scalar t0,
    scalar dt,
    scalarField& values
)
{
    (void)t0;
    (void)dt;

    const fvMesh& mesh = domain.mesh();

    if (!reportedConductivitySource_)
    {
        if (Pstream::master())
        {
            Info<< "pseudoECG: using conductivity field '"
                << domain.conductivity().name() << "' on mesh '" << mesh.name()
                << "'" << nl << endl;
        }

        reportedConductivitySource_ = true;
    }

    if (!leadVectorsCalculated_)
    {
        calculateLeadVectors(domain);
    }

    const tmp<volVectorField> tgradVm = fvc::grad(domain.Vm());
    const vectorField& gradVm = tgradVm().primitiveField();

    const label nElectrodes = domain.electrodePositions().size();

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

        values[electrodeI] = ecgVal;
    }

    for (label electrodeI = 0; electrodeI < nElectrodes; ++electrodeI)
    {
        reduce(values[electrodeI], sumOp<scalar>());
    }

    if (sigmaE_ > VSMALL)
    {
        const scalar norm =
            1.0 / (4.0 * constant::mathematical::pi * sigmaE_);
        forAll(values, eI)
        {
            values[eI] *= norm;
        }
    }

    if (!hasOwnSampling_)
    {
        return;
    }

    const scalar currentTime = mesh.time().value();

    domain.recordVerification(currentTime, values);

    if (currentTime + SMALL < startTime_ || currentTime > endTime_ + SMALL)
    {
        return;
    }

    if (currentTime + SMALL < nextSampleTime_)
    {
        return;
    }

    if (!outputPtr_.valid())
    {
        outputPtr_ =
            ecgModelIO::openTimeSeries
            (
                mesh.time().globalPath()/"postProcessing",
                "pseudoECG.dat",
                domain.electrodeNames()
            );
    }

    ecgModelIO::writeRow(outputPtr_.ref(), currentTime, values);
    nextSampleTime_ += deltaT_;
}

} // End namespace Foam

// ************************************************************************* //
