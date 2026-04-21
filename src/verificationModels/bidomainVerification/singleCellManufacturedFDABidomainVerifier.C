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

#include "singleCellManufacturedFDABidomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "ionicModel.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "addToRunTimeSelectionTable.H"

#include <tuple>

namespace Foam
{

defineTypeNameAndDebug(singleCellManufacturedFDABidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    singleCellManufacturedFDABidomainVerifier,
    dictionary
);

namespace
{

label requireFieldIndex(const wordList& names, const word& name, const char* phase)
{
    forAll(names, i)
    {
        if (names[i] == name)
        {
            return i;
        }
    }

    FatalErrorInFunction
        << "Required " << phase << " field '" << name
        << "' is missing from configured hook fields " << names
        << exit(FatalError);

    return -1;
}


bool shouldReportManufacturedErrors(const volScalarField& Vm)
{
    const Time& time = Vm.mesh().time();
    const scalar t = time.value();
    const scalar dt = time.deltaTValue();
    const scalar endTime = time.endTime().value();

    return t + 0.5*dt >= endTime;
}


label globalManufacturedCellCount(const fvMesh& mesh)
{
    label totalCells = mesh.nCells();
    reduce(totalCells, sumOp<label>());
    return totalCells;
}


label structuredCellsPerDirection(const label totalCells, const label dimension)
{
    if (totalCells <= 0 || dimension <= 0)
    {
        return 0;
    }

    if (dimension == 1)
    {
        return totalCells;
    }

    return max
    (
        label(1),
        label(Foam::pow(scalar(totalCells), 1.0/scalar(dimension)) + 0.5)
    );
}


scalar structuredManufacturedDx(const label nPerDirection)
{
    if (nPerDirection <= 0)
    {
        return 0.0;
    }

    return 1.0/scalar(nPerDirection);
}


void computeAndWriteManufacturedErrors
(
    const scalarField& Vm,
    const scalarField& u1m,
    const scalarField& u2m,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const scalar dimension,
    const int N,
    const scalar dx,
    const scalar dt,
    const int nsteps,
    const bool useExplicitAlgorithm,
    const fileName& outFilename
)
{
    scalarField Vex;
    scalarField u1ex;
    scalarField u2ex;
    scalarField u3ex;
    computeManufacturedV(Vex, x, y, z, t, dimension);
    computeManufacturedU(u1ex, u2ex, u3ex, x, y, z, t, dimension);

    auto computeNorms = [&](const scalarField& num, const scalarField& exact)
    {
        scalar sumAbs = 0.0;
        scalar sumSq = 0.0;
        scalar maxAbs = 0.0;

        forAll(num, i)
        {
            const scalar diff = Foam::mag(num[i] - exact[i]);
            sumAbs += diff;
            sumSq += diff*diff;
            if (diff > maxAbs)
            {
                maxAbs = diff;
            }
        }

        reduce(maxAbs, maxOp<scalar>());
        reduce(sumAbs, sumOp<scalar>());
        reduce(sumSq, sumOp<scalar>());

        label n = num.size();
        reduce(n, sumOp<label>());

        return std::tuple<scalar, scalar, scalar>
        (
            sumAbs/scalar(n),
            Foam::sqrt(sumSq/scalar(n)),
            maxAbs
        );
    };

    const auto [L1_V, L2_V, Linf_V] = computeNorms(Vm, Vex);
    const auto [L1_u1, L2_u1, Linf_u1] = computeNorms(u1m, u1ex);
    const auto [L1_u2, L2_u2, Linf_u2] = computeNorms(u2m, u2ex);

    if (!Pstream::master())
    {
        return;
    }

    Info << "\nSingle-cell bidomain manufactured-solution error summary (t = "
         << t << "):" << nl
         << "-------------------------------------------------\n"
         << "Field     L1-error       L2-error       Linf-error\n"
         << "Vm     " << L1_V << "   " << L2_V << "   " << Linf_V << nl
         << "u1     " << L1_u1 << "   " << L2_u1 << "   " << Linf_u1 << nl
         << "u2     " << L1_u2 << "   " << L2_u2 << "   " << Linf_u2 << nl
         << "-------------------------------------------------\n" << endl;

    OFstream out(outFilename);
    out << "Single-cell bidomain manufactured-solution error summary (t = "
        << t << "):\n";
    out << "Field     L1-error       L2-error       Linf-error\n";
    out << "Vm     " << L1_V << "   " << L2_V << "   " << Linf_V << "\n";
    out << "u1     " << L1_u1 << "   " << L2_u1 << "   " << Linf_u1 << "\n";
    out << "u2     " << L1_u2 << "   " << L2_u2 << "   " << Linf_u2 << "\n";
    out << "-------------------------------------------------\n\n";

    out << "\nSimulation summary:\n";
    out << "-------------------\n";
    out << "Number of cells (N)   = " << N << "\n";
    out << "Solver type           = "
        << (useExplicitAlgorithm ? "Explicit" : "Implicit") << "\n";
    out << "Grid spacing (dx)     = " << dx << "\n";
    out << "Time step (dt)        = " << dt << "\n";
    out << "Number of steps       = " << nsteps << "\n";
    out << "Final simulation time = " << t << "\n";
    out << "-------------------\n\n";
}

} // End anonymous namespace


singleCellManufacturedFDABidomainVerifier::
singleCellManufacturedFDABidomainVerifier
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

wordList singleCellManufacturedFDABidomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList singleCellManufacturedFDABidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2"});
}


bool singleCellManufacturedFDABidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void singleCellManufacturedFDABidomainVerifier::preProcess
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
            << "singleCellManufacturedFDABidomainVerifier preProcess expected "
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


void singleCellManufacturedFDABidomainVerifier::postProcess
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
            << "singleCellManufacturedFDABidomainVerifier postProcess expected "
            << requiredNames.size() << " scratch fields " << requiredNames
            << " but received " << fields.size()
            << exit(FatalError);
    }

    const label u1Field =
        requireFieldIndex(requiredNames, "u1", "postProcess requirement");
    const label u2Field =
        requireFieldIndex(requiredNames, "u2", "postProcess requirement");

    const fvMesh& mesh = Vm.mesh();
    const vectorField& centres = mesh.C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));

    const scalarField& VmValues = Vm.primitiveField();
    const scalarField& u1Values = fields[u1Field].primitiveField();
    const scalarField& u2Values = fields[u2Field].primitiveField();
    const scalar t = mesh.time().value();
    const scalar dt = mesh.time().deltaTValue();
    const label dimension = model.geometricDimension();
    const label totalCells = globalManufacturedCellCount(mesh);
    const label nPerDirection =
        structuredCellsPerDirection(totalCells, dimension);
    const scalar dx = structuredManufacturedDx(nPerDirection);
    const label nSteps = max(label(0), mesh.time().timeIndex());
    const fileName outputDir(mesh.time().path()/"postProcessing");
    mkDir(outputDir);
    const word dimName =
        dimension == 1 ? "1D"
      : dimension == 2 ? "2D"
      : dimension == 3 ? "3D"
                       : "unknown";
    const fileName outputFile
    (
        outputDir
      / (
            dimName
          + "_"
          + Foam::name(nPerDirection)
          + "_cells_"
          + word(useExplicitAlgorithm_ ? "explicit" : "implicit")
          + ".dat"
        )
    );

    computeAndWriteManufacturedErrors
    (
        VmValues,
        u1Values,
        u2Values,
        X,
        Y,
        Z,
        t,
        dimension,
        nPerDirection,
        dx,
        dt,
        nSteps,
        useExplicitAlgorithm_,
        outputFile
    );

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
