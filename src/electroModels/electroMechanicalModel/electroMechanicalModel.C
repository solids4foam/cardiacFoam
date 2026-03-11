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
#include "ionicModelIO.H"
#include "OSspecific.H"
#include "stimulusIO.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace electroModels {
defineTypeNameAndDebug(electroMechanicalModel, 0);
addToRunTimeSelectionTable(physicsModel, electroMechanicalModel, physicsModel);
} // namespace electroModels
} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::physicsModel>
Foam::electroModels::electroMechanicalModel::New(Time &runTime,
                                                 const word &region) {
  return autoPtr<physicsModel>(
      new electroMechanicalModel(typeName, runTime, region));
}

Foam::electroModels::electroMechanicalModel::electroMechanicalModel(
    Time &runTime, const word &region)
    : electroMechanicalModel(typeName, runTime, region) {}

Foam::electroModels::electroMechanicalModel::electroMechanicalModel(
    const word &type, Time &runTime, const word &region)
    : physicsModel(type, runTime),
      electroModelPtr_(electroModel::New(runTime, region)),
      activeTensionModelPtr_(), TaPtr_(), lambdaPtr_(),
      activeTensionOutFields_(), activeTensionOutputPtr_() {
  const dictionary &props = electroModelPtr_->electroProperties();

  if (props.found("activeTensionModel")) {
    const dictionary &atSubDict = props.subDict("activeTensionModel");
    const word atModelType(atSubDict.lookup("activeTensionModel"));

    // Merge parent dict so ODE solver keys (solver, initialODEStep, maxSteps,
    // etc.) are available in the subdict without duplicating them in the case.
    dictionary atDict(atSubDict);
    atDict.merge(props);

    // Keep the model selector as a primitive entry.
    // merge(props) injects the parent "activeTensionModel" sub-dictionary,
    // which would otherwise replace this key with a dictionary entry.
    atDict.add("activeTensionModel", atModelType, true);

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

    // Optional exported active-tension fields for volumetric output
    const wordList atNames = activeTensionModelPtr_->exportedFieldNames();
    activeTensionOutFields_.setSize(atNames.size());
    forAll(atNames, i) {
      activeTensionOutFields_.set(
          i, new volScalarField(IOobject(atNames[i], runTime.timeName(),
                                         electroModelPtr_->mesh(),
                                         IOobject::NO_READ, IOobject::AUTO_WRITE),
                                electroModelPtr_->mesh(), dimless,
                                "zeroGradient"));
    }

    // Single-cell trace output for active-tension model variables
    if (electroModelPtr_->mesh().nCells() == 1) {
      const fileName outputDir(runTime.path() / "postProcessing");
      mkDir(outputDir);

      const word ionicModelType(props.lookup("ionicModel"));
      const word tissueName(props.lookup("tissue"));
      const word stimSuffix(stimulusIO::protocolSuffix(props));
      const fileName outFile =
          outputDir / (ionicModelType + "_" + tissueName + "_" + stimSuffix +
                       "_activeTension_" + atModelType + ".txt");
      activeTensionOutputPtr_.reset(new OFstream(outFile));

      OFstream &output = activeTensionOutputPtr_.ref();
      output.setf(std::ios::fixed);
      output.precision(7);
      activeTensionModelPtr_->writeHeader(output);
    }
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

    if (activeTensionOutputPtr_.valid()) {
      const scalar t1 = runTime().value();
      const scalar t0 = t1 - runTime().deltaTValue();
      if (ionicModelIO::shouldWriteStep(t0, t1, electroModelPtr_->electroProperties(),
                                        false)) {
        activeTensionModelPtr_->write(t1, activeTensionOutputPtr_.ref());
      }
    }
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
  if (activeTensionOutFields_.size()) {
    activeTensionModelPtr_->exportStates(activeTensionOutFields_);
  }
  if (TaPtr_.valid())
    TaPtr_->write();
}

void Foam::electroModels::electroMechanicalModel::end() {
  electroModelPtr_->end();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
