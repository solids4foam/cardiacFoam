/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pseudoECGElectro.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "OSspecific.H"
#include "ecgModelIO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pseudoECGElectro, 0);
addToRunTimeSelectionTable(ecgModel, pseudoECGElectro, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

pseudoECGElectro::pseudoECGElectro(const volScalarField &Vm,
                                   const dictionary &dict,
                                   const volTensorField *conductivityPtr)
    : ecgModel(Vm, dict, conductivityPtr), outputPtr_() {
  const fileName outDir(mesh_.time().path() / "postProcessing");
  outputPtr_ =
      ecgModelIO::openTimeSeries(outDir, "pseudoECG.dat", electrodeNames_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pseudoECGElectro::compute() {
  // Gima-Rudy dipole:
  //   phi_pseudo(P) = -sum_c [ (Gi . grad(Vm))_c . r_vec * V_c / |r|^3 ]
  //   r_vec = C_c - P  (no sigmaT factor, consistent with original impl)

  const tmp<volVectorField> tgradVm = fvc::grad(Vm_);
  const vectorField &gradVm = tgradVm().primitiveField();

  const scalarField &Vols = mesh_.V();
  const vectorField &Ctrs = mesh_.C().primitiveField();
  const tensorField &GiI = Gi().primitiveField();

  const label nE = electrodePositions_.size();

  List<scalar> localSums(nE, scalar(0));

  forAll(Ctrs, cI) {
    const vector dipole = (GiI[cI] & gradVm[cI]) * Vols[cI];

    for (label pI = 0; pI < nE; pI++) {
      const vector r_vec = Ctrs[cI] - electrodePositions_[pI];
      const scalar r = mag(r_vec);
      if (r > VSMALL) {
        localSums[pI] += (dipole & r_vec) / (r * r * r);
      }
    }
  }

  // Parallel reduction (no sign flip: consistent with original accumulation)
  for (label pI = 0; pI < nE; pI++) {
    reduce(localSums[pI], sumOp<scalar>());
  }

  // Write to pseudoECG.dat at outputTime
  if (mesh_.time().outputTime()) {
    ecgModelIO::writeRow(outputPtr_.ref(), mesh_.time().value(), localSums);
  }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
