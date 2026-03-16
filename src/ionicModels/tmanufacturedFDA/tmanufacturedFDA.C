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

#include <math.h>
#include "tmanufacturedFDA.H"
#include "tmanufacturedFDA_2014.H"
#include "tmanufacturedFDA_2014Names.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "ionicSelector.H"

#include "tmanufacturedFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
defineTypeNameAndDebug(tmanufacturedFDA, 0);
addToRunTimeSelectionTable(ionicModel, tmanufacturedFDA, dictionary);
} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tmanufacturedFDA::tmanufacturedFDA(const dictionary &dict,
                                         const label num,
                                         const scalar initialDeltaT,
                                         const Switch solveVmWithinODESolver)
    : ionicModel(dict, num, initialDeltaT, solveVmWithinODESolver),
      STATES_(num), STATES_OLD_(num), CONSTANTS_(NUM_CONSTANTS, 0.0),
      ALGEBRAIC_(num), RATES_(num) {
  // Dimension read from dict: "dimension 2D;" — avoids mesh dependency in base
  const word dimWord(dict.lookup("dimension"));
  setTissue(ionicSelector::dimensionFlag(dimWord));

  forAll(STATES_, integrationPtI) {
    STATES_.set(integrationPtI, new scalarField(NUM_STATES, 0.0));
    STATES_OLD_.set(integrationPtI, new scalarField(NUM_STATES, 0.0));
    ALGEBRAIC_.set(integrationPtI, new scalarField(NUM_ALGEBRAIC, 0.0));
    RATES_.set(integrationPtI, new scalarField(NUM_STATES, 0.0));

    // Initialise the constants (repeatedly! it's ok...) and the rates and
    // states
    tmanufacturedFDAinitConsts(CONSTANTS_.data(), RATES_[integrationPtI].data(),
                               STATES_[integrationPtI].data(), tissue());
    STATES_OLD_[integrationPtI] = STATES_[integrationPtI];
  }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tmanufacturedFDA::~tmanufacturedFDA() {}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::tmanufacturedFDA::supportedDimensions() const

{
  return {"1D", "2D", "3D"};
}

void Foam::tmanufacturedFDA::preProcess(volScalarField &Vm,
                                        PtrList<volScalarField> &outFields) {
  if (outFields.size() < 3) {
    FatalErrorInFunction
        << "tmanufacturedFDA requires at least 3 outFields (u1, u2, u3)"
        << exit(FatalError);
  }

  volScalarField &u1m = outFields[0];
  volScalarField &u2m = outFields[1];
  volScalarField &u3m = outFields[2];

  const volVectorField &C = Vm.mesh().C();
  scalarField X = C.component(vector::X);
  scalarField Y = C.component(vector::Y);
  scalarField Z = C.component(vector::Z);
  const scalar t = 0.0;

  computeManufacturedV(Vm, X, Y, Z, t, tissue());
  computeManufacturedU(u1m, u2m, u3m, X, Y, Z, t, tissue());

  forAll(STATES_, i) {
    scalarField &S = STATES_[i];
    S[V] = Vm[i];
    S[u1] = u1m[i];
    S[u2] = u2m[i];
    S[u3] = u3m[i];

    STATES_OLD_[i] = S;
  }

  Vm.correctBoundaryConditions();
  u1m.correctBoundaryConditions();
  u2m.correctBoundaryConditions();
  u3m.correctBoundaryConditions();

  Info << "Boundary conditions corrected." << endl;
}

void Foam::tmanufacturedFDA::postProcess(
    const volScalarField &Vm, const PtrList<volScalarField> &outFields) {
  // MS logic for error computation could go here if needed.
  // Currently, error reporting is handled via external call or handled
  // by the msHandler in v1. In v2, we should move it here.
  // However, I don't see computeAndPrintErrors called in this file yet.
  // I'll check tmanufacturedFields.H for it.
}

void Foam::tmanufacturedFDA::solveODE(const scalar stepStartTime,
                                      const scalar deltaT,
                                      const scalarField &Vm, scalarField &Im) {
  // STATES_ = STATES_OLD_;
  const scalar tStart = stepStartTime;
  const scalar tEnd = tStart + deltaT;
  const label monitorCell = 0;

  forAll(STATES_, integrationPtI) {
    scalarField &S = STATES_[integrationPtI];
    scalarField &A = ALGEBRAIC_[integrationPtI];
    scalarField &R = RATES_[integrationPtI];

    scalar &h = ionicModel::step()[integrationPtI];

    S[V] = Vm[integrationPtI];
    h = min(h, deltaT);

    if (integrationPtI == monitorCell) {
      debugPrintFields(integrationPtI, tStart, tEnd, h);
    }

    odeSolver().solve(tStart, tEnd, S, h);

    ::tmanufacturedFDAcomputeVariables(tEnd, CONSTANTS_.data(), R.data(),
                                       S.data(), A.data(), tissue(),
                                       solveVmWithinODESolver());
    if (integrationPtI == monitorCell) {
      debugPrintFields(integrationPtI, tStart, tEnd, h);
    }

    Im[integrationPtI] = A[Iion];
  }
}
void Foam::tmanufacturedFDA::derivatives(const scalar t, const scalarField &y,
                                         scalarField &dydt) const {
  scalarField ALG(NUM_ALGEBRAIC, 0.0);

  ::tmanufacturedFDAcomputeVariables(
      t, CONSTANTS_.data(), dydt.data(), const_cast<scalarField &>(y).data(),
      ALG.data(), tissue(), solveVmWithinODESolver());
}

void Foam::tmanufacturedFDA::updateStatesOld() const {
  saveStateSnapshot(STATES_, STATES_OLD_);
}

void Foam::tmanufacturedFDA::resetStatesToStatesOld() const {
  restoreStateSnapshot(STATES_, STATES_OLD_);
}

// ------------------------------------------------------------------------- //
//  Writing logic for 1D-3D manufactured solution

// exporting always 3 fields for manufactured computation
Foam::wordList Foam::tmanufacturedFDA::exportedFieldNames() const {

  Foam::wordList names;
  names.append("u1");
  names.append("u2");
  names.append("u3");
  return names;
}

// ************************************************************************* //
