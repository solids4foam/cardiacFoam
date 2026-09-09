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

#include "manufacturedAnisotropicMonodomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "boundBox.H"
#include "ionicModel.H"
#include "monodomainVerification/manufacturedAnisotropicMonodomainReference.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedAnisotropicMonodomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedAnisotropicMonodomainVerifier,
    dictionary
);


manufacturedAnisotropicMonodomainVerifier::
manufacturedAnisotropicMonodomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    errorsReported_(false),
    beta_(0.0),
    betaInitialised_(false),
    conductivityValidated_(false)
{}


void manufacturedAnisotropicMonodomainVerifier::initialiseBeta
(
    const ionicModel& model
)
{
    const wordList constantNames = model.constantVariableNames();
    const scalarField constantValues = model.constantVariableValues();

    if (constantNames.size() != constantValues.size())
    {
        FatalErrorInFunction
            << "Ionic model " << model.type()
            << " exposes inconsistent constant metadata: "
            << constantNames.size() << " names but "
            << constantValues.size() << " values."
            << exit(FatalError);
    }

    label betaIndex = -1;
    forAll(constantNames, i)
    {
        if (constantNames[i] == "Beta")
        {
            betaIndex = i;
            break;
        }
    }

    if (betaIndex < 0)
    {
        FatalErrorInFunction
            << "manufacturedAnisotropicMonodomainVerifier requires the "
            << "manufactured ionic constant 'Beta', but ionic model "
            << model.type() << " exposes constants " << constantNames << "."
            << exit(FatalError);
    }

    beta_ = constantValues[betaIndex];
    betaInitialised_ = true;
}


void manufacturedAnisotropicMonodomainVerifier::validateConfiguration
(
    const ionicModel& model,
    const fvMesh& mesh
) const
{
    if
    (
        model.verificationFamily() != "monodomainFDAManufactured"
     || model.geometricDimension() != 3
    )
    {
        FatalErrorInFunction
            << "manufacturedAnisotropicMonodomainVerifier requires the 3D "
            << "monodomainFDAManufactured ionic model. Received ionic model "
            << model.type() << " with verification family '"
            << model.verificationFamily() << "' and geometric dimension "
            << model.geometricDimension() << "."
            << exit(FatalError);
    }

    const boundBox bounds(mesh.points(), true);
    const vector unitMin(vector::zero);
    const vector unitMax(vector::one);
    const scalar cubeTolerance = 1e-10;

    if
    (
        mag(bounds.min() - unitMin) > cubeTolerance
     || mag(bounds.max() - unitMax) > cubeTolerance
    )
    {
        FatalErrorInFunction
            << "manufacturedAnisotropicMonodomainVerifier is defined on "
            << "the unit cube [0,1]^3 so that its homogeneous-flux boundary "
            << "condition is exact. Mesh bounds are " << bounds << "."
            << exit(FatalError);
    }
}


void manufacturedAnisotropicMonodomainVerifier::validateConductivity
(
    const volTensorField& conductivity
) const
{
    const tensorField& sigma = conductivity.primitiveField();

    tensor sigmaSum(tensor::zero);
    forAll(sigma, cellI)
    {
        sigmaSum += sigma[cellI];
    }
    reduce(sigmaSum, sumOp<tensor>());

    label nCells = sigma.size();
    reduce(nCells, sumOp<label>());

    if (nCells <= 0)
    {
        FatalErrorInFunction
            << "Cannot validate anisotropic conductivity on an empty mesh."
            << exit(FatalError);
    }

    const tensor meanSigma = sigmaSum/scalar(nCells);
    scalar maxDeviation = 0.0;
    forAll(sigma, cellI)
    {
        maxDeviation = max(maxDeviation, mag(sigma[cellI] - meanSigma));
    }
    reduce(maxDeviation, maxOp<scalar>());

    const scalar scale = max(scalar(1), mag(meanSigma));
    const scalar tolerance = 1e-10*scale;

    if (maxDeviation > tolerance)
    {
        FatalErrorInFunction
            << "manufacturedAnisotropicMonodomainVerifier assumes a "
            << "spatially constant conductivity tensor. Maximum cellwise "
            << "deviation from the global mean is " << maxDeviation
            << ", exceeding tolerance " << tolerance << "."
            << exit(FatalError);
    }

    const scalar symmetryError = mag(meanSigma - meanSigma.T());
    if (symmetryError > tolerance)
    {
        FatalErrorInFunction
            << "The manufactured anisotropic conductivity must be symmetric. "
            << "The global mean tensor is " << meanSigma
            << " and ||K-K^T|| = " << symmetryError << "."
            << exit(FatalError);
    }

    const tensor symmetricSigma = 0.5*(meanSigma + meanSigma.T());
    const scalar leadingMinor1 = symmetricSigma.xx();
    const scalar leadingMinor2 =
        symmetricSigma.xx()*symmetricSigma.yy()
      - Foam::sqr(symmetricSigma.xy());
    const scalar leadingMinor3 = det(symmetricSigma);

    if
    (
        leadingMinor1 <= tolerance
     || leadingMinor2 <= Foam::sqr(tolerance)
     || leadingMinor3 <= Foam::pow(tolerance, 3)
    )
    {
        FatalErrorInFunction
            << "The manufactured anisotropic conductivity must be symmetric "
            << "positive definite. Sylvester principal minors are ("
            << leadingMinor1 << ' ' << leadingMinor2 << ' '
            << leadingMinor3 << ") for tensor " << symmetricSigma << "."
            << exit(FatalError);
    }

    const scalar maxOffDiagonal = max
    (
        mag(symmetricSigma.xy()),
        max(mag(symmetricSigma.xz()), mag(symmetricSigma.yz()))
    );

    if (maxOffDiagonal <= tolerance)
    {
        FatalErrorInFunction
            << "The manufactured anisotropic conductivity must contain a "
            << "genuine rotation (at least one non-zero off-diagonal "
            << "component). Received " << symmetricSigma << "."
            << exit(FatalError);
    }
}


wordList manufacturedAnisotropicMonodomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList
manufacturedAnisotropicMonodomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2"});
}


bool manufacturedAnisotropicMonodomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedAnisotropicMonodomainVerifier::preProcess
(
    ionicModel& model,
    volScalarField& Vm,
    PtrList<volScalarField>& fields
)
{
    const wordList names = preProcessFieldNames(model);

    if (fields.size() != names.size())
    {
        FatalErrorInFunction
            << "manufacturedAnisotropicMonodomainVerifier preProcess "
            << "expected " << names.size() << " preProcess fields " << names
            << " but received " << fields.size()
            << exit(FatalError);
    }

    validateConfiguration(model, Vm.mesh());
    initialiseBeta(model);

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    volScalarField& u1m = fields[u1Field];
    volScalarField& u2m = fields[u2Field];
    volScalarField& u3m = fields[u3Field];

    const vectorField& centres = Vm.mesh().C().primitiveField();
    const scalar t = Vm.mesh().time().value();

    scalarField& VmI = Vm.primitiveFieldRef();
    scalarField& u1I = u1m.primitiveFieldRef();
    scalarField& u2I = u2m.primitiveFieldRef();
    scalarField& u3I = u3m.primitiveFieldRef();

    computeAnisotropicManufacturedV(VmI, centres, t);
    computeAnisotropicManufacturedU(u1I, u2I, u3I, centres, t);

    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    model.importFields(Vm, names, fields);
}


void manufacturedAnisotropicMonodomainVerifier::addManufacturedPdeSource
(
    volScalarField& sourceField,
    const volTensorField& conductivity,
    const scalar evaluationTime
) const
{
    if (!betaInitialised_)
    {
        FatalErrorInFunction
            << "The manufactured source was requested before ionic constant "
            << "Beta was initialised by preProcess."
            << exit(FatalError);
    }

    if (&sourceField.mesh() != &conductivity.mesh())
    {
        FatalErrorInFunction
            << "Source field and conductivity field belong to different "
            << "finite-volume meshes."
            << exit(FatalError);
    }

    if (!conductivityValidated_)
    {
        validateConductivity(conductivity);
        conductivityValidated_ = true;
    }

    scalarField manufacturedSource;
    computeAnisotropicManufacturedPdeSource
    (
        manufacturedSource,
        sourceField.mesh().C().primitiveField(),
        conductivity.primitiveField(),
        evaluationTime,
        beta_
    );

    sourceField.primitiveFieldRef() += manufacturedSource;
    sourceField.correctBoundaryConditions();
}


void manufacturedAnisotropicMonodomainVerifier::postProcess
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

    const wordList requiredNames = requiredPostProcessFieldNames(model);

    if (fields.size() != requiredNames.size())
    {
        FatalErrorInFunction
            << "manufacturedAnisotropicMonodomainVerifier postProcess "
            << "expected " << requiredNames.size() << " scratch fields "
            << requiredNames << " but received " << fields.size()
            << exit(FatalError);
    }

    const label u1Idx = requireFieldIndex(requiredNames, "u1", "postProcess");
    const label u2Idx = requireFieldIndex(requiredNames, "u2", "postProcess");

    const fvMesh& mesh = Vm.mesh();
    const vectorField& centres = mesh.C().primitiveField();
    const scalar t = mesh.time().value();
    const scalar dt = mesh.time().deltaTValue();
    const label totalCells = globalManufacturedCellCount(mesh);
    const label nPerDirection = structuredCellsPerDirection(totalCells, 3);
    const scalar dx = structuredManufacturedDx(nPerDirection);
    const label nSteps = max(label(0), mesh.time().timeIndex());

    scalarField VmExact, u1Exact, u2Exact, u3Exact;
    computeAnisotropicManufacturedV(VmExact, centres, t);
    computeAnisotropicManufacturedU
    (
        u1Exact,
        u2Exact,
        u3Exact,
        centres,
        t
    );

    const auto VmNorms = computeNorms(mesh, Vm.primitiveField(), VmExact);
    const auto u1Norms =
        computeNorms(mesh, fields[u1Idx].primitiveField(), u1Exact);
    const auto u2Norms =
        computeNorms(mesh, fields[u2Idx].primitiveField(), u2Exact);

    const fileName outputDir(mesh.time().globalPath()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile
    (
        outputDir
      / (
            word("3D_") + Foam::name(nPerDirection) + "_cells.dat"
        )
    );

    if (Pstream::master())
    {
        Info << "\nRotated-anisotropy manufactured-solution error summary "
             << "(t = " << t << "):" << nl
             << "-------------------------------------------------" << nl
             << "Field     L1-error       L2-error       Linf-error" << nl
             << "Vm     " << VmNorms.first().first() << "   "
             << VmNorms.first().second() << "   " << VmNorms.second() << nl
             << "u1     " << u1Norms.first().first() << "   "
             << u1Norms.first().second() << "   " << u1Norms.second() << nl
             << "u2     " << u2Norms.first().first() << "   "
             << u2Norms.first().second() << "   " << u2Norms.second() << nl
             << "-------------------------------------------------" << endl;

        OFstream out(outputFile);
        out << "Rotated-anisotropy manufactured-solution error summary "
            << "(t = " << t << "):\n";
        out << "Field     L1-error       L2-error       Linf-error\n";
        out << "Vm     " << VmNorms.first().first() << "   "
            << VmNorms.first().second() << "   " << VmNorms.second() << "\n";
        out << "u1     " << u1Norms.first().first() << "   "
            << u1Norms.first().second() << "   " << u1Norms.second() << "\n";
        out << "u2     " << u2Norms.first().first() << "   "
            << u2Norms.first().second() << "   " << u2Norms.second() << "\n";
        out << "-------------------------------------------------\n\n";

        out << "Simulation summary:\n";
        out << "-------------------\n";
        out << "Number of cells (N)   = " << nPerDirection << "\n";

        out << "Grid spacing (dx)     = " << dx << "\n";
        out << "Time step (dt)        = " << dt << "\n";
        out << "Number of steps       = " << nSteps << "\n";
        out << "Final simulation time = " << t << "\n";
        out << "Beta                  = " << beta_ << "\n";
        out << "-------------------\n\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
