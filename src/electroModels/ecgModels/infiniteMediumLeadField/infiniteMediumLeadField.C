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

#include "infiniteMediumLeadField.H"

#include "ecgDomain.H"
#include "PstreamReduceOps.H"
#include "mathematicalConstants.H"

namespace Foam
{

infiniteMediumLeadField::infiniteMediumLeadField(const dictionary& dict)
:
    sigmaExtracellular_(0.0),
    leadVectorsCalculated_(false),
    leadVectors_()
{
    if (!dict.found("sigmaExtracellular"))
    {
        FatalIOErrorInFunction(dict)
            << "ECG domain '" << dict.dictName() << "' has no "
            << "'sigmaExtracellular': the scalar extracellular conductivity "
            << "of the infinite medium is required."
            << exit(FatalIOError);
    }

    sigmaExtracellular_ = dict.get<scalar>("sigmaExtracellular");

    if (sigmaExtracellular_ <= 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "ECG domain '" << dict.dictName() << "': "
            << "'sigmaExtracellular' must be positive, got "
            << sigmaExtracellular_ << "."
            << exit(FatalIOError);
    }
}


scalar infiniteMediumLeadField::scale() const
{
    return 1.0 / (4.0 * constant::mathematical::pi * sigmaExtracellular_);
}


void infiniteMediumLeadField::calculateLeadVectors(const ecgDomain& domain)
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


void infiniteMediumLeadField::evaluate
(
    const ecgDomain& domain,
    const vectorField& gradVm,
    scalarField& values
)
{
    if (!leadVectorsCalculated_)
    {
        calculateLeadVectors(domain);
    }

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

        reduce(ecgVal, sumOp<scalar>());

        values[electrodeI] = ecgVal;
    }

    const scalar norm = scale();
    forAll(values, eI)
    {
        values[eI] *= norm;
    }
}

} // End namespace Foam

// ************************************************************************* //
