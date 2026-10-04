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

#include "manufacturedConormalMonodomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "boundBox.H"
#include "ionicModel.H"
#include "mathematicalConstants.H"
#include "monodomainVerification/manufacturedConormalMonodomainReference.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedConormalMonodomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedConormalMonodomainVerifier,
    dictionary
);


manufacturedConormalMonodomainVerifier::
manufacturedConormalMonodomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    errorsReported_(false),
    beta_(0.0),
    betaInitialised_(false),
    profileInitialised_(false),
    profile_{complex(0, 0), complex(0, 0), complex(0, 0), 0, 0}
{}


void manufacturedConormalMonodomainVerifier::initialiseBeta
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
            << "manufacturedConormalMonodomainVerifier requires the "
            << "manufactured ionic constant 'Beta', but ionic model "
            << model.type() << " exposes constants " << constantNames << "."
            << exit(FatalError);
    }

    beta_ = constantValues[betaIndex];
    betaInitialised_ = true;
}


void manufacturedConormalMonodomainVerifier::validateConfiguration
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
            << "manufacturedConormalMonodomainVerifier requires the 3D "
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
            << "manufacturedConormalMonodomainVerifier is defined on "
            << "the unit cube [0,1]^3. Mesh bounds are " << bounds << "."
            << exit(FatalError);
    }
}


void manufacturedConormalMonodomainVerifier::initialiseProfile
(
    const volTensorField& conductivity
) const
{
    if (profileInitialised_)
    {
        return;
    }

    const tensor meanSigma =
        uniformTensor(conductivity.primitiveField(), typeName);

    profile_ = conormalWallProfile
    (
        meanSigma,
        2*constant::mathematical::pi,
        2*constant::mathematical::pi
    );
    profileInitialised_ = true;
}


wordList manufacturedConormalMonodomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList
manufacturedConormalMonodomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2"});
}


bool manufacturedConormalMonodomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedConormalMonodomainVerifier::preProcess
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
            << "manufacturedConormalMonodomainVerifier preProcess "
            << "expected " << names.size() << " preProcess fields " << names
            << " but received " << fields.size()
            << exit(FatalError);
    }

    validateConfiguration(model, Vm.mesh());
    initialiseBeta(model);

    const volTensorField& conductivity =
        Vm.mesh().lookupObject<volTensorField>("conductivity");
    initialiseProfile(conductivity);

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

    computeConormalManufacturedV(VmI, centres, t, profile_);
    computeConormalManufacturedU(u1I, u2I, u3I, centres, t, profile_);

    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    model.importFields(Vm, names, fields);
}


void manufacturedConormalMonodomainVerifier::addManufacturedPdeSource
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

    initialiseProfile(conductivity);

    scalarField manufacturedSource;
    computeConormalManufacturedPdeSource
    (
        manufacturedSource,
        sourceField.mesh().C().primitiveField(),
        conductivity.primitiveField(),
        evaluationTime,
        beta_,
        profile_
    );

    sourceField.primitiveFieldRef() += manufacturedSource;
    sourceField.correctBoundaryConditions();
}


void manufacturedConormalMonodomainVerifier::postProcess
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
            << "manufacturedConormalMonodomainVerifier postProcess "
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
    computeConormalManufacturedV(VmExact, centres, t, profile_);
    computeConormalManufacturedU
    (
        u1Exact,
        u2Exact,
        u3Exact,
        centres,
        t,
        profile_
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
        Info << "\nConormal manufactured-solution error summary "
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
        out << "Conormal manufactured-solution error summary "
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
