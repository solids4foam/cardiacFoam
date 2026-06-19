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

#include "coupled1D3DMonodomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "electroDomainCouplingEndpoints.H"
#include "electroDomainInterface.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(coupled1D3DMonodomainVerifier, 0);
addToRunTimeSelectionTable
(
    couplingVerificationModel,
    coupled1D3DMonodomainVerifier,
    dictionary
);


coupled1D3DMonodomainVerifier::coupled1D3DMonodomainVerifier
(
    const dictionary& dict
)
:
    couplingVerificationModel(dict),
    diagnosticsWritten_(false)
{}


void coupled1D3DMonodomainVerifier::preProcess
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain
)
{
    (void)primaryDomain;
    (void)secondaryDomain;
}


void coupled1D3DMonodomainVerifier::postProcess
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain
)
{
    if (diagnosticsWritten_)
    {
        return;
    }

    diagnosticsWritten_ = true;

    const fvMesh& mesh = primaryDomain.mesh();
    const volScalarField& sourceField = primaryDomain.sourceField();
    const scalarField& sourceValues = sourceField.primitiveField();
    const scalarField& volumes = mesh.V();

    scalar sourceMin = GREAT;
    scalar sourceMax = -GREAT;
    scalar sumAbs = 0.0;
    scalar sumSq = 0.0;
    scalar totalCurrent = 0.0;
    scalar totalAbsCurrent = 0.0;
    label nonZeroCells = 0;

    forAll(sourceValues, i)
    {
        const scalar source = sourceValues[i];
        sourceMin = min(sourceMin, source);
        sourceMax = max(sourceMax, source);
        sumAbs += mag(source);
        sumSq += source*source;
        totalCurrent += source*volumes[i];
        totalAbsCurrent += mag(source)*volumes[i];

        if (mag(source) > SMALL)
        {
            ++nonZeroCells;
        }
    }

    reduce(sourceMin, minOp<scalar>());
    reduce(sourceMax, maxOp<scalar>());
    reduce(sumAbs, sumOp<scalar>());
    reduce(sumSq, sumOp<scalar>());
    reduce(totalCurrent, sumOp<scalar>());
    reduce(totalAbsCurrent, sumOp<scalar>());
    reduce(nonZeroCells, sumOp<label>());

    const label totalCells = globalManufacturedCellCount(mesh);
    const scalar denom = max(scalar(1), scalar(totalCells));
    const scalar sourceL1 = sumAbs/denom;
    const scalar sourceL2 = Foam::sqrt(sumSq/denom);

    if (Pstream::master())
    {
        mkDir(mesh.time().globalPath()/"verification");
        autoPtr<OFstream> osPtr
        (
            new OFstream
            (
                mesh.time().globalPath()/"verification"
              / "coupled1D3DMonodomain_diagnostics.csv"
            )
        );

        if (osPtr.valid())
        {
            *osPtr
                << "t,sourceMin,sourceMax,sourceL1,sourceL2,"
                << "totalCurrent,totalAbsCurrent,nonZeroCells,totalCells" << nl;

            *osPtr
                << mesh.time().value() << ','
                << sourceMin << ','
                << sourceMax << ','
                << sourceL1 << ','
                << sourceL2 << ','
                << totalCurrent << ','
                << totalAbsCurrent << ','
                << nonZeroCells << ','
                << totalCells << nl;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
