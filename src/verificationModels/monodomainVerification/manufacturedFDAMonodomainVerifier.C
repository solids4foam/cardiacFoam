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

#include "manufacturedFDAMonodomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "ionicModel.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedFDAMonodomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedFDAMonodomainVerifier,
    dictionary
);

manufacturedFDAMonodomainVerifier::manufacturedFDAMonodomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    useExplicitAlgorithm_
    (
        dict.lookupOrDefault<word>("solutionAlgorithm", "implicit")
     == "explicit"
    ),
    errorsReported_(false)
{}

wordList manufacturedFDAMonodomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList manufacturedFDAMonodomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2"});
}


bool manufacturedFDAMonodomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedFDAMonodomainVerifier::preProcess
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
            << "manufacturedFDAMonodomainVerifier preProcess expected "
            << names.size() << " preProcess fields " << names
            << " but received " << fields.size()
            << exit(FatalError);
    }

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    volScalarField& u1m = fields[u1Field];
    volScalarField& u2m = fields[u2Field];
    volScalarField& u3m = fields[u3Field];

    const vectorField& centres = Vm.mesh().C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));
    const scalar t = Vm.mesh().time().value();

    scalarField& VmI = Vm.primitiveFieldRef();
    scalarField& u1I = u1m.primitiveFieldRef();
    scalarField& u2I = u2m.primitiveFieldRef();
    scalarField& u3I = u3m.primitiveFieldRef();

    const label dimension = model.geometricDimension();

    computeManufacturedV(VmI, X, Y, Z, t, dimension);
    computeManufacturedU(u1I, u2I, u3I, X, Y, Z, t, dimension);

    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    model.importFields(Vm, names, fields);
}


void manufacturedFDAMonodomainVerifier::postProcess
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
            << "manufacturedFDAMonodomainVerifier postProcess expected "
            << requiredNames.size() << " scratch fields " << requiredNames
            << " but received " << fields.size()
            << exit(FatalError);
    }

    const label u1Idx = requireFieldIndex(requiredNames, "u1", "postProcess");
    const label u2Idx = requireFieldIndex(requiredNames, "u2", "postProcess");

    const fvMesh& mesh = Vm.mesh();
    const vectorField& centres = mesh.C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));

    const scalar t = mesh.time().value();
    const scalar dt = mesh.time().deltaTValue();
    const label dimension = model.geometricDimension();
    const label totalCells = globalManufacturedCellCount(mesh);
    const label nPerDirection = structuredCellsPerDirection(totalCells, dimension);
    const scalar dx = structuredManufacturedDx(nPerDirection);
    const label nSteps = max(label(0), mesh.time().timeIndex());

    scalarField VmExact, u1Exact, u2Exact, u3Exact;
    computeManufacturedV(VmExact, X, Y, Z, t, dimension);
    computeManufacturedU(u1Exact, u2Exact, u3Exact, X, Y, Z, t, dimension);

    const auto VmNorms = computeNorms(Vm.primitiveField(), VmExact);
    const auto u1Norms = computeNorms(fields[u1Idx].primitiveField(), u1Exact);
    const auto u2Norms = computeNorms(fields[u2Idx].primitiveField(), u2Exact);

    const fileName outputDir(mesh.time().path()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile
    (
        outputDir
      / (
            dimensionName(dimension)
          + "_"
          + Foam::name(nPerDirection)
          + "_cells_"
          + word(useExplicitAlgorithm_ ? "explicit" : "implicit")
          + ".dat"
        )
    );

    if (Pstream::master())
    {
        Info << "\nSimulation summary:\n"
             << "-------------------\n"
             << "Number of cells (N)   = " << nPerDirection << nl
             << "Solver type           = "
             << (useExplicitAlgorithm_ ? "Explicit" : "Implicit") << nl
             << "Grid spacing (dx)     = " << dx << nl
             << "Time step (dt)        = " << dt << nl
             << "Number of steps       = " << nSteps << nl
             << "Final simulation time = " << t << nl
             << "-------------------\n";

        Info << "\nManufactured-solution error summary (t = " << t << "):" << nl
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
        out << "Manufactured-solution error summary (t = " << t << "):\n";
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
        out << "Solver type           = "
            << (useExplicitAlgorithm_ ? "Explicit" : "Implicit") << "\n";
        out << "Grid spacing (dx)     = " << dx << "\n";
        out << "Time step (dt)        = " << dt << "\n";
        out << "Number of steps       = " << nSteps << "\n";
        out << "Final simulation time = " << t << "\n";
        out << "-------------------\n\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
