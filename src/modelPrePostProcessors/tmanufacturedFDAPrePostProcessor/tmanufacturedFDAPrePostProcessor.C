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

#include "tmanufacturedFDAPrePostProcessor.H"

#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "OSspecific.H"
#include "tmanufacturedExactFields.H"
#include <tuple>

namespace Foam {
namespace electroModels {

namespace {

wordList readConfiguredNames(const dictionary &dict, const word &subDictName,
                             const word &entryName,
                             const wordList &defaults) {
  if (!dict.found(subDictName)) {
    return defaults;
  }

  const dictionary &hookDict = dict.subDict(subDictName);
  return hookDict.lookupOrDefault<wordList>(entryName, defaults);
}

label requireFieldIndex(const wordList &names, const word &name,
                        const char *phase) {
  forAll(names, i) {
    if (names[i] == name) {
      return i;
    }
  }

  FatalErrorInFunction << "Required " << phase << " field '" << name
                       << "' is missing from configured hook fields " << names
                       << exit(FatalError);
  return -1;
}

bool shouldReportManufacturedErrors(const volScalarField &Vm) {
  const Time &time = Vm.mesh().time();
  const scalar t = time.value();
  const scalar dt = time.deltaTValue();
  const scalar endTime = time.endTime().value();

  return t + 0.5 * dt >= endTime;
}

label globalManufacturedCellCount(const fvMesh &mesh) {
  label totalCells = mesh.nCells();
  reduce(totalCells, sumOp<label>());
  return totalCells;
}

label structuredCellsPerDirection(const label totalCells,
                                  const label dimension) {
  if (totalCells <= 0 || dimension <= 0) {
    return 0;
  }

  if (dimension == 1) {
    return totalCells;
  }

  return max(
      label(1),
      label(Foam::pow(scalar(totalCells), 1.0 / scalar(dimension)) + 0.5));
}

scalar structuredManufacturedDx(const label nPerDirection) {
  if (nPerDirection <= 0) {
    return 0.0;
  }

  // Manufactured tutorial meshes are structured unit domains in the active
  // directions, so dx is fully determined by the global logical resolution.
  return 1.0 / scalar(nPerDirection);
}

void computeAndWriteManufacturedErrors(
    const scalarField &Vm, const scalarField &u1m, const scalarField &u2m,
    const scalarField &x, const scalarField &y, const scalarField &z,
    const scalar t, const scalar dimension, const int N, const scalar dx,
    const scalar dt, const int nsteps, const bool useExplicitAlgorithm,
    const fileName &outFilename) {
  scalarField Vex, u1ex, u2ex, u3ex;
  computeManufacturedV(Vex, x, y, z, t, dimension);
  computeManufacturedU(u1ex, u2ex, u3ex, x, y, z, t, dimension);

  auto computeNorms = [&](const scalarField &num, const scalarField &exact) {
    scalar sumAbs = 0.0, sumSq = 0.0, maxAbs = 0.0;
    forAll(num, i) {
      const scalar diff = Foam::mag(num[i] - exact[i]);
      sumAbs += diff;
      sumSq += diff * diff;
      if (diff > maxAbs) {
        maxAbs = diff;
      }
    }

    reduce(maxAbs, maxOp<scalar>());
    reduce(sumAbs, sumOp<scalar>());
    reduce(sumSq, sumOp<scalar>());
    label n = num.size();
    reduce(n, sumOp<int>());

    const scalar L1 = sumAbs / n;
    const scalar L2 = Foam::sqrt(sumSq / n);
    const scalar Linf = maxAbs;
    return std::tuple<scalar, scalar, scalar>(L1, L2, Linf);
  };

  const auto [L1_V, L2_V, Linf_V] = computeNorms(Vm, Vex);
  const auto [L1_u1, L2_u1, Linf_u1] = computeNorms(u1m, u1ex);
  const auto [L1_u2, L2_u2, Linf_u2] = computeNorms(u2m, u2ex);

  if (!Pstream::master()) {
    return;
  }

  Info << "\nSimulation summary:\n";
  Info << "-------------------\n";
  Info << "Number of cells (N)   = " << N << nl;
  Info << "Solver type           = "
       << (useExplicitAlgorithm ? "Explicit" : "Implicit") << nl;
  Info << "Grid spacing (dx)     = " << dx << nl;
  Info << "Time step (dt)        = " << dt << nl;
  Info << "Number of steps       = " << nsteps << nl;
  Info << "Final simulation time = " << t << nl;
  Info << "-------------------\n";
  Info << "\nManufactured-solution error summary (t = " << t << "):" << nl
       << "-------------------------------------------------\n"
       << "Field     L1-error       L2-error       Linf-error\n"
       << "Vm     " << L1_V << "   " << L2_V << "   " << Linf_V << nl
       << "u1     " << L1_u1 << "   " << L2_u1 << "   " << Linf_u1 << nl
       << "u2     " << L1_u2 << "   " << L2_u2 << "   " << Linf_u2 << nl
       << "-------------------------------------------------\n" << endl;

  OFstream fout(outFilename);
  fout << "Manufactured-solution error summary (t = " << t << "):\n";
  fout << "Field     L1-error       L2-error       Linf-error\n";
  fout << "Vm     " << L1_V << "   " << L2_V << "   " << Linf_V << "\n";
  fout << "u1     " << L1_u1 << "   " << L2_u1 << "   " << Linf_u1 << "\n";
  fout << "u2     " << L1_u2 << "   " << L2_u2 << "   " << Linf_u2 << "\n";
  fout << "-------------------------------------------------\n\n";

  fout << "\nSimulation summary:\n";
  fout << "-------------------\n";
  fout << "Number of cells (N)   = " << N << "\n";
  fout << "Solver type           = "
       << (useExplicitAlgorithm ? "Explicit" : "Implicit") << "\n";
  fout << "Grid spacing (dx)     = " << dx << "\n";
  fout << "Time step (dt)        = " << dt << "\n";
  fout << "Number of steps       = " << nsteps << "\n";
  fout << "Final simulation time = " << t << "\n";
  fout << "-------------------\n\n";
}

} // End anonymous namespace

defineTypeNameAndDebug(tmanufacturedFDAPrePostProcessor, 0);
addToRunTimeSelectionTable(electroModelsPrePostProcessor,
                           tmanufacturedFDAPrePostProcessor, dictionary);

tmanufacturedFDAPrePostProcessor::tmanufacturedFDAPrePostProcessor(
    const dictionary &dict)
    : electroModelsPrePostProcessor(dict),
      useExplicitAlgorithm_(
          dict.lookupOrDefault<word>("solutionAlgorithm", "implicit") ==
          "explicit"),
      errorsReported_(false) {}

wordList
tmanufacturedFDAPrePostProcessor::preProcessFieldNames(
    const ionicModel &) const {
  return readConfiguredNames(dict(), "solverHookFields", "preProcess",
                             wordList({"u1", "u2", "u3"}));
}

wordList
tmanufacturedFDAPrePostProcessor::requiredPostProcessFieldNames(
    const ionicModel &) const {
  return wordList({"u1", "u2"});
}

bool tmanufacturedFDAPrePostProcessor::shouldPostProcess(
    const ionicModel &, const volScalarField &Vm) const {
  return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}

void tmanufacturedFDAPrePostProcessor::preProcess(
    ionicModel &model, volScalarField &Vm, PtrList<volScalarField> &fields) {
  const wordList names = preProcessFieldNames(model);

  if (fields.size() != names.size()) {
    FatalErrorInFunction
        << "tmanufacturedFDA preProcess expected " << names.size()
        << " preProcess fields " << names << " but received " << fields.size()
        << exit(FatalError);
  }

  const label u1Field = requireFieldIndex(names, "u1", "preProcess");
  const label u2Field = requireFieldIndex(names, "u2", "preProcess");
  const label u3Field = requireFieldIndex(names, "u3", "preProcess");

  volScalarField &u1m = fields[u1Field];
  volScalarField &u2m = fields[u2Field];
  volScalarField &u3m = fields[u3Field];

  const volVectorField &C = Vm.mesh().C();
  const vectorField &centres = C.primitiveField();
  scalarField X(centres.component(vector::X));
  scalarField Y(centres.component(vector::Y));
  scalarField Z(centres.component(vector::Z));
  const scalar t = Vm.mesh().time().value();

  scalarField &VmI = Vm.primitiveFieldRef();
  scalarField &u1I = u1m.primitiveFieldRef();
  scalarField &u2I = u2m.primitiveFieldRef();
  scalarField &u3I = u3m.primitiveFieldRef();

  const label dimension = model.tissue();

  computeManufacturedV(VmI, X, Y, Z, t, dimension);
  computeManufacturedU(u1I, u2I, u3I, X, Y, Z, t, dimension);

  Vm.correctBoundaryConditions();
  u1m.correctBoundaryConditions();
  u2m.correctBoundaryConditions();
  u3m.correctBoundaryConditions();

  model.importFields(Vm, names, fields);
}

void tmanufacturedFDAPrePostProcessor::postProcess(
    const ionicModel &model, const volScalarField &Vm,
    const PtrList<volScalarField> &fields) {
  if (!shouldPostProcess(model, Vm)) {
    return;
  }

  const wordList names = model.exportedFieldNames();
  const wordList requiredNames = requiredPostProcessFieldNames(model);

  if (fields.size() != names.size()) {
    FatalErrorInFunction
        << "tmanufacturedFDA postProcess expected " << names.size()
        << " exported fields " << names << " but received " << fields.size()
        << exit(FatalError);
  }

  const label u1Field =
      requireFieldIndex(requiredNames, "u1", "postProcess requirement");
  const label u2Field =
      requireFieldIndex(requiredNames, "u2", "postProcess requirement");
  const label exportedU1Field =
      requireFieldIndex(names, requiredNames[u1Field], "postProcess export");
  const label exportedU2Field =
      requireFieldIndex(names, requiredNames[u2Field], "postProcess export");

  const fvMesh &mesh = Vm.mesh();
  const volVectorField &C = mesh.C();
  const vectorField &centres = C.primitiveField();
  scalarField X(centres.component(vector::X));
  scalarField Y(centres.component(vector::Y));
  scalarField Z(centres.component(vector::Z));

  const scalarField &VmValues = Vm.primitiveField();
  const scalarField &u1Values = fields[exportedU1Field].primitiveField();
  const scalarField &u2Values = fields[exportedU2Field].primitiveField();
  const scalar t = mesh.time().value();
  const scalar dt = mesh.time().deltaTValue();
  const label dimension = model.tissue();
  const label totalCells = globalManufacturedCellCount(mesh);
  const label nPerDirection =
      structuredCellsPerDirection(totalCells, dimension);
  const scalar dx = structuredManufacturedDx(nPerDirection);
  const label nSteps = max(label(0), mesh.time().timeIndex());
  const fileName outputDir(mesh.time().path() / "postProcessing");
  mkDir(outputDir);
  const word dimName =
      (dimension == 1)   ? "1D"
      : (dimension == 2) ? "2D"
      : (dimension == 3) ? "3D"
                         : "unknown";
  const fileName outputFile(
      outputDir / (dimName + "_" + Foam::name(nPerDirection) + "_cells_" +
                   word(useExplicitAlgorithm_ ? "explicit" : "implicit") +
                   ".dat"));

  computeAndWriteManufacturedErrors(
      VmValues, u1Values, u2Values, X, Y, Z, t, dimension, nPerDirection, dx,
      dt, nSteps, useExplicitAlgorithm_, outputFile);

  errorsReported_ = true;
}

} // End namespace electroModels
} // End namespace Foam
