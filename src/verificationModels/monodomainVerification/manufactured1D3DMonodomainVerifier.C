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

#include "manufactured1D3DMonodomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "ionicModel.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufactured1D3DMonodomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufactured1D3DMonodomainVerifier,
    dictionary
);


manufactured1D3DMonodomainVerifier::manufactured1D3DMonodomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    sourceFieldPtr_(nullptr),
    reportBaseline3DReference_
    (
        dict.lookupOrDefault<Switch>("reportBaseline3DReference", false)
    ),
    diagnosticsWritten_(false)
{}


void manufactured1D3DMonodomainVerifier::bindSourceField
(
    const volScalarField& sourceField
)
{
    sourceFieldPtr_ = &sourceField;
}


wordList manufactured1D3DMonodomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList();
}


wordList manufactured1D3DMonodomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList();
}


bool manufactured1D3DMonodomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !diagnosticsWritten_ && shouldReportManufacturedErrors(Vm);
}


void manufactured1D3DMonodomainVerifier::preProcess
(
    ionicModel&,
    volScalarField&,
    PtrList<volScalarField>& fields
)
{
    if (!fields.empty())
    {
        FatalErrorInFunction
            << "manufactured1D3DMonodomainVerifier does not expect "
            << "preProcess fields but received " << fields.size()
            << exit(FatalError);
    }
}


void manufactured1D3DMonodomainVerifier::postProcess
(
    const ionicModel& model,
    const volScalarField& Vm,
    const PtrList<volScalarField>& fields
)
{
    if (!shouldPostProcess(model, Vm))
    {
        return;
    }

    if (!fields.empty())
    {
        FatalErrorInFunction
            << "manufactured1D3DMonodomainVerifier does not expect "
            << "postProcess fields but received " << fields.size()
            << exit(FatalError);
    }

    if (!sourceFieldPtr_)
    {
        FatalErrorInFunction
            << "manufactured1D3DMonodomainVerifier requires a bound "
            << "myocardium sourceField."
            << exit(FatalError);
    }

    const fvMesh& mesh = Vm.mesh();
    const scalarField& sourceValues = sourceFieldPtr_->primitiveField();
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
    const label dimension = model.geometricDimension();
    const label nPerDirection = structuredCellsPerDirection(totalCells, dimension);
    const scalar dx = structuredManufacturedDx(nPerDirection);
    const scalar denom = max(scalar(1), scalar(totalCells));
    const scalar sourceL1 = sumAbs/denom;
    const scalar sourceL2 = Foam::sqrt(sumSq/denom);
    const scalar sourceLinf = max(mag(sourceMin), mag(sourceMax));
    const scalar t = mesh.time().value();
    const scalar dt = mesh.time().deltaTValue();

    scalar baselineVmL1 = 0.0;
    scalar baselineVmL2 = 0.0;
    scalar baselineVmLinf = 0.0;

    if (reportBaseline3DReference_)
    {
        const vectorField& centres = mesh.C().primitiveField();
        scalarField X(centres.component(vector::X));
        scalarField Y(centres.component(vector::Y));
        scalarField Z(centres.component(vector::Z));
        scalarField VmExact;
        computeManufacturedV(VmExact, X, Y, Z, t, dimension);

        const auto VmNorms = computeNorms(Vm.primitiveField(), VmExact);
        baselineVmL1 = VmNorms.first().first();
        baselineVmL2 = VmNorms.first().second();
        baselineVmLinf = VmNorms.second();
    }

    const fileName outputDir(mesh.time().path()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile
    (
        outputDir
      / (
            "manufactured1D3D_"
          + dimensionName(dimension)
          + "_"
          + Foam::name(nPerDirection)
          + "_cells.dat"
        )
    );

    if (Pstream::master())
    {
        Info<< "\n1D-3D manufactured coupling diagnostics (t = " << t << "):"
            << nl
            << "-------------------------------------------------" << nl
            << "dt                       = " << dt << nl
            << "dx                       = " << dx << nl
            << "active source cells      = " << nonZeroCells << nl
            << "source[min,max]          = [" << sourceMin << ", "
            << sourceMax << "]" << nl
            << "source L1/L2/Linf        = "
            << sourceL1 << ", " << sourceL2 << ", " << sourceLinf << nl
            << "integrated source        = " << totalCurrent << nl
            << "integrated |source|      = " << totalAbsCurrent << nl;

        if (reportBaseline3DReference_)
        {
            Info<< "baseline 3D Vm L1/L2/Linf = "
                << baselineVmL1 << ", " << baselineVmL2 << ", "
                << baselineVmLinf << nl;
        }

        Info<< "-------------------------------------------------" << endl;

        OFstream out(outputFile);
        out << "1D-3D manufactured coupling diagnostics (t = " << t << ")\n";
        out << "dt " << dt << "\n";
        out << "dx " << dx << "\n";
        out << "activeSourceCells " << nonZeroCells << "\n";
        out << "sourceMin " << sourceMin << "\n";
        out << "sourceMax " << sourceMax << "\n";
        out << "sourceL1 " << sourceL1 << "\n";
        out << "sourceL2 " << sourceL2 << "\n";
        out << "sourceLinf " << sourceLinf << "\n";
        out << "integratedSource " << totalCurrent << "\n";
        out << "integratedAbsSource " << totalAbsCurrent << "\n";

        if (reportBaseline3DReference_)
        {
            out << "baseline3DVmL1 " << baselineVmL1 << "\n";
            out << "baseline3DVmL2 " << baselineVmL2 << "\n";
            out << "baseline3DVmLinf " << baselineVmLinf << "\n";
        }
    }

    diagnosticsWritten_ = true;
}

} // End namespace Foam

// ************************************************************************* //
