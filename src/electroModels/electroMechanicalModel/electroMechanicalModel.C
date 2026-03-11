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

#include "electroMechanicalModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace electroModels {
defineTypeNameAndDebug(electroMechanicalModel, 0);
addToRunTimeSelectionTable(physicsModel, electroMechanicalModel, physicsModel);
} // namespace electroModels
} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electroModels::electroMechanicalModel::electroMechanicalModel(
    const word &type, Time &runTime, const word &region)
    : physicsModel(type, runTime),
      electroModelPtr_(electroModel::New(runTime, region)),
      activeTensionModelPtr_(), TaPtr_(), lambdaPtr_() {
  const dictionary &props = electroModelPtr_->electroProperties();

  if (props.found("activeTensionModel")) {
    // Merge parent dict so ODE solver keys (solver, initialODEStep, maxSteps,
    // etc.) are available in the subdict without duplicating them in the case.
    dictionary atDict(props.subDict("activeTensionModel"));
    atDict.merge(props);
    activeTensionModelPtr_ =
        activeTensionModel::New(atDict, electroModelPtr_->mesh().nCells());

    // Connect signal provider (ionic model) to active tension
    const CouplingSignalProvider *provider = electroModelPtr_->provider();
    if (provider) {
      activeTensionModelPtr_->setCouplingSignalProvider(*provider);
    }

    // Active tension field
    TaPtr_.reset(new volScalarField(
        IOobject("Ta", runTime.timeName(), electroModelPtr_->mesh(),
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        electroModelPtr_->mesh(), dimensionedScalar("Ta", dimPressure, 0.0),
        "zeroGradient"));

    // Lambda field (stretch) - typically read from solid simulation
    lambdaPtr_.reset(new volScalarField(
        IOobject("lambda", runTime.timeName(), electroModelPtr_->mesh(),
                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        electroModelPtr_->mesh(), dimensionedScalar("lambda", dimless, 1.0),
        "zeroGradient"));
  }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroModels::electroMechanicalModel::~electroMechanicalModel() {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electroModels::electroMechanicalModel::setDeltaT(Time &runTime) {
  electroModelPtr_->setDeltaT(runTime);
}

bool Foam::electroModels::electroMechanicalModel::evolve() {
  bool status = electroModelPtr_->evolve();

  if (activeTensionModelPtr_.valid()) {
    activeTensionModelPtr_->calculateTension(
        runTime().value(), runTime().deltaTValue(), lambdaPtr_->internalField(),
        TaPtr_->primitiveFieldRef());
    TaPtr_->correctBoundaryConditions();
  }

  return status;
}

bool Foam::electroModels::electroMechanicalModel::read() {
  bool status = electroModelPtr_->read();

  if (activeTensionModelPtr_.valid()) {
    // Re-read active tension if needed
  }

  return status;
}

void Foam::electroModels::electroMechanicalModel::writeFields(
    const Time &runTime) {
  electroModelPtr_->writeFields(runTime);
  if (TaPtr_.valid())
    TaPtr_->write();
}

void Foam::electroModels::electroMechanicalModel::end() {
  electroModelPtr_->end();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
