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

    const auto VmNorms = computeNorms(mesh, Vm.primitiveField(), VmExact);
    const auto u1Norms =
        computeNorms(mesh, fields[u1Idx].primitiveField(), u1Exact);
    const auto u2Norms =
        computeNorms(mesh, fields[u2Idx].primitiveField(), u2Exact);

    // Bulk and boundary-adjacent parts of the Vm error and of its cell gradient
    const scalar S = Foam::sqrt(1.0 + t);
    vectorField gradVmExact(centres.size());
    forAll(centres, cellI)
    {
        const point& c = centres[cellI];
        gradVmExact[cellI] =
            S*pseudoECGManufacturedGradF(c.x(), c.y(), c.z(), dimension);
    }
    const tmp<volVectorField> tgradVm(fvc::grad(Vm));
    const tmp<volVectorField> tgradRecon
    (
        exactFieldGradient
        (
            mesh,
            Vm.name(),
            [&](const point& p) { return S*manufacturedFDAFAt(p, dimension); }
        )
    );
    const auto splitRows = bulkBoundaryRows
    (
        mesh, "Vm", Vm.primitiveField(), VmExact,
        tgradVm().primitiveField(), tgradRecon().primitiveField(), gradVmExact
    );

    const fileName outputDir(mesh.time().globalPath()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile
    (
        outputDir
      / (
            dimensionName(dimension) + "_" + Foam::name(nPerDirection) + "_cells.dat"
        )
    );

    if (Pstream::master())
    {
        errorTable table;
        table.titlePrefix = "Manufactured-solution error summary (t = ";
        table.t = t;
        table.addRow("Vm     ", VmNorms);
        table.addRow("u1     ", u1Norms);
        table.addRow("u2     ", u2Norms);
        table.splitRows = splitRows;

        Info << "\nSimulation summary:\n"
             << "-------------------\n";
        writeStructuredRunLines(Info, nPerDirection, dx, dt, nSteps, t);
        Info << "-------------------\n";

        writeErrorTableConsole(Info, table);

        OFstream out(outputFile);
        writeErrorTableFile(out, table);
        out << "-------------------------------------------------\n\n";

        out << "Simulation summary:\n";
        out << "-------------------\n";
        writeStructuredRunLines(out, nPerDirection, dx, dt, nSteps, t);
        out << "-------------------\n\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
