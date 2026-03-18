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
      STATES_(num), CONSTANTS_(NUM_CONSTANTS, 0.0), ALGEBRAIC_(num),
      RATES_(num) {
  setTissue(ionicSelector::selectDimension(dict, supportedDimensions()));

  forAll(STATES_, integrationPtI) {
    STATES_.set(integrationPtI, new scalarField(NUM_STATES, 0.0));
    ALGEBRAIC_.set(integrationPtI, new scalarField(NUM_ALGEBRAIC, 0.0));
    RATES_.set(integrationPtI, new scalarField(NUM_STATES, 0.0));

    // Initialise the constants (repeatedly! it's ok...) and the rates and
    // states
    tmanufacturedFDAinitConsts(CONSTANTS_.data(), RATES_[integrationPtI].data(),
                               STATES_[integrationPtI].data(), tissue());
  }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tmanufacturedFDA::~tmanufacturedFDA() {}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::tmanufacturedFDA::supportedDimensions() const

{
  return {"1D", "2D", "3D"};
}

void Foam::tmanufacturedFDA::solveODE(const scalar stepStartTime,
                                      const scalar deltaT,
                                      const scalarField &Vm, scalarField &Im) {
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

    // The generic monodomain solver expects ionic current normalized by Cm.
    Im[integrationPtI] = A[Iion] / CONSTANTS_[Cm];
  }
}
void Foam::tmanufacturedFDA::derivatives(const scalar t, const scalarField &y,
                                         scalarField &dydt) const {
  scalarField ALG(NUM_ALGEBRAIC, 0.0);

  ::tmanufacturedFDAcomputeVariables(
      t, CONSTANTS_.data(), dydt.data(), const_cast<scalarField &>(y).data(),
      ALG.data(), tissue(), solveVmWithinODESolver());
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
