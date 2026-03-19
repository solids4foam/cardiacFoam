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

#include "monoDomainElectro.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "IOmanip.H"
#include "stimulusIO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

namespace electroModels {

namespace {

bool hasCustomPrePostProcessor(const electroModelsPrePostProcessor &processor) {
  return processor.type() != "electroModelsPrePostProcessor";
}

} // End anonymous namespace

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(monoDomainElectro, 0);
addToRunTimeSelectionTable(electroModel, monoDomainElectro, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

tmp<volTensorField> monoDomainElectro::initialiseConductivityTensor(
    const word &fieldName, const dictionary &dict, const fvMesh &mesh) {
  tmp<volTensorField> tresult(new volTensorField(
      IOobject(fieldName, mesh.time().timeName(), mesh,
               IOobject::READ_IF_PRESENT, IOobject::NO_WRITE),
      mesh,
      dimensionedTensor("zero",
                        pow3(dimTime) * sqr(dimCurrent) / (dimMass * dimVolume),
                        tensor::zero)));
  volTensorField &result = tresult.ref();

  if (!result.headerOk()) {
    Info << "\n"
         << fieldName << " not found on disk, using value from " << dict.name()
         << nl << endl;

    result = dimensionedTensor(
        dimensionedSymmTensor(
            fieldName, pow3(dimTime) * sqr(dimCurrent) / (dimMass * dimVolume),
            dict) &
        tensor(I));

    if (result.size() > 0) {
      Info << fieldName << " tensor (cell 0): " << result[0] << nl;
    }
  } else {
    Info << fieldName << " field read from " << mesh.time().timeName() << nl
         << endl;
  }

  return tresult;
}

tmp<volTensorField> monoDomainElectro::initialiseConductivity() const {
  return initialiseConductivityTensor("conductivity", cardiacProperties_,
                                      mesh());
}

void monoDomainElectro::readVm() { Vm_.read(); }

void monoDomainElectro::updateExternalStimulusCurrent(
    volScalarField &externalStimulusCurrent,
    const List<stimulus> &externalStimulus, const scalar t0) const {
  scalarField &externalStimulusCurrentI = externalStimulusCurrent;
  externalStimulusCurrentI = 0.0;

  forAll(externalStimulus, stimI) {
    const stimulus &curStim = externalStimulus[stimI];

    const scalar tStart = curStim.startTime_;
    if (t0 < tStart || t0 > (tStart + curStim.duration_)) {
      continue;
    }

    const labelList &stimulusCellIDs = curStim.cellIDs_;
    forAll(stimulusCellIDs, cI) {
      const label id = stimulusCellIDs[cI];
      externalStimulusCurrentI[id] += curStim.intensity_;
    }
  }

  externalStimulusCurrent.correctBoundaryConditions();
}

void monoDomainElectro::updateActivationTime(volScalarField &activationTime,
                                             boolList &calculateActivationTime,
                                             const volScalarField &Vm) const {
  // Update activationTime field
  const scalarField &VmI = Vm;
  const scalarField &VmOldI = Vm.oldTime();
  boolList &calculateActivationTimeI = calculateActivationTime;
  scalarField &activationTimeI = activationTime.primitiveFieldRef();
  const scalar oldTime = runTime().value() - runTime().deltaTValue();
  const scalar deltaT = runTime().deltaTValue();
  forAll(activationTimeI, cellI) {
    if (calculateActivationTimeI[cellI]) {
      if (VmI[cellI] > SMALL) {
        calculateActivationTimeI[cellI] = false;

        // Linearly interpolate for more accuracy
        const scalar w = (0.0 - VmOldI[cellI]) / (VmI[cellI] - VmOldI[cellI]);

        activationTimeI[cellI] = oldTime + w * deltaT;
      }
    }
  }

  activationTime.correctBoundaryConditions();
}

bool monoDomainElectro::evolveExplicit() {
  if (time().timeIndex() == 1) {
    Info << "Solving the mono-domain equation for Vm using an explicit "
         << "approach" << nl << setw(20) << "Simulation Time" << setw(20)
         << "Clock Time" << setw(20) << "Min Vm" << setw(20) << "Max Vm"
         << endl;
  }

  physicsModel::printInfo() =
      bool(time().timeIndex() % infoFrequency_ == 0 ||
           mag(time().value() - time().endTime().value()) < SMALL);

  if (physicsModel::printInfo()) {
    Info << setw(20) << time().value() << setw(20) << time().elapsedClockTime()
         << setw(20) << min(Vm_).value() << setw(20) << max(Vm_).value()
         << endl;

    physicsModel::printInfo() = false;
  }

  // Disable OpenFOAM linear solver output
  const label debugOrg = SolverPerformance<scalar>::debug;
  SolverPerformance<scalar>::debug = 0;

  // Current time step
  const scalar dt = runTime().deltaTValue();

  // Old time
  const scalar t0 = runTime().value() - dt;

  // 1) Set the external stimulus
  updateExternalStimulusCurrent(externalStimulusCurrent_, externalStimulus_,
                                t0);

  // 2) Advance ionic model in time using the current Vm state
  ionicModelPtr_->solveODE(t0, dt,
                           Vm_,
                           Iion_);
  Iion_.correctBoundaryConditions();

  // 3) Solve the Vm equation using the updated ionic current
  solve(chi_ * Cm_ * fvm::ddt(Vm_) == fvc::laplacian(conductivity_, Vm_) -
                                          chi_ * Cm_ * Iion_ +
                                          externalStimulusCurrent_);

  // 5) Update post-processing fields

  updateActivationTime(activationTime_, calculateActivationTime_, Vm_);

  // Re-enable OpenFOAM linear solver output
  SolverPerformance<scalar>::debug = debugOrg;

  return true;
}

bool monoDomainElectro::evolveImplicit() {
  if (time().timeIndex() == 1) {
    Info << "Solving the mono-domain equation for Vm using an implicit "
         << "approach" << endl;
  }

  // Current time step
  const scalar dt = runTime().deltaTValue();

  // Old time
  const scalar t0 = runTime().value() - dt;

  // 1) Set the external stimulus
  updateExternalStimulusCurrent(externalStimulusCurrent_, externalStimulus_,
                                t0);

  // 2) Advance ionic model in time using the current Vm state
  ionicModelPtr_->solveODE(t0, dt,
                           Vm_,
                           Iion_);
  Iion_.correctBoundaryConditions();

  // 3) Solve the Vm equation implicitly using the updated ionic current
  while (pimple().loop()) {
    solve(chi_ * Cm_ * fvm::ddt(Vm_) == fvm::laplacian(conductivity_, Vm_) -
                                            chi_ * Cm_ * Iion_ +
                                            externalStimulusCurrent_);
  }

  // 5) Update post-processing fields

  updateActivationTime(activationTime_, calculateActivationTime_, Vm_);

  return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

monoDomainElectro::monoDomainElectro(Time &runTime, const word &region)
    : monoDomainElectro(typeName, runTime, region) {}

monoDomainElectro::monoDomainElectro(const word &type, Time &runTime,
                                     const word &region)
    : electroModel(type, runTime, region),
      Vm_(IOobject("Vm", runTime.timeName(), mesh(), IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE),
          mesh(), dimensionedScalar("Vm", dimVoltage, -0.084), "zeroGradient"),
      Iion_(IOobject("ionicCurrent", runTime.timeName(), mesh(),
                     IOobject::NO_READ, IOobject::AUTO_WRITE),
            mesh(), dimensionedScalar("zero", dimVoltage / dimTime, 0.0),
            "zeroGradient"),
      externalStimulusCurrent_(
          IOobject("externalStimulusCurrent", runTime.timeName(), mesh(),
                   IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
          mesh(), dimensionedScalar("zero", dimCurrent / dimVolume, 0.0),
          "zeroGradient"),
      externalStimulus_(),
      cardiacProperties_(electroProperties()),
      chi_("chi", dimArea / dimVolume, cardiacProperties_),
      Cm_("Cm", dimCurrent * dimTime / (dimVoltage * dimArea),
          cardiacProperties_),
      conductivity_(initialiseConductivity()),
      activationTime_(IOobject("activationTime", runTime.timeName(), mesh(),
                               IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                      mesh(), dimensionedScalar("zero", dimTime, 0.0),
                      "zeroGradient"),
      calculateActivationTime_(mesh().nCells(), true), outFields_(),
      ionicModelPtr_(ionicModel::New(ionicProperties(),
                                     mesh().nCells(),
                                     runTime.deltaTValue())),
      prePostProcessor_(electroModelsPrePostProcessor::New(
          ionicModelPtr_->type(), ionicProperties())),
      preProcessFieldNames_(prePostProcessor_->preProcessFieldNames(*ionicModelPtr_)),
      preProcessFields_(),
      setDeltaT_(true),
      infoFrequency_(
          electroProperties().lookupOrDefault<label>("infoFrequency", 1)) {
  const Switch allowIonicStimulusInMonodomain =
      electroProperties().lookupOrDefault<Switch>(
          "allowIonicStimulusInMonodomain", false);

  const StimulusProtocol &ionicStim = ionicModelPtr_->stimulusProtocol();
  if (stimulusIO::hasActiveStimulus(ionicStim)) {
    if (!allowIonicStimulusInMonodomain) {
      FatalErrorInFunction
          << "Detected active ionic-model stimulus protocol while using "
          << "monoDomainElectro." << nl
          << "This can unintentionally combine ionic and PDE external "
          << "stimulation." << nl << "Either disable ionic protocol keys "
          << "(stim_*, nstim*) or set "
          << "allowIonicStimulusInMonodomain true to override."
          << abort(FatalError);
    }

    Info << "Warning: active ionic-model stimulus is enabled in "
         << "monoDomainElectro (allowIonicStimulusInMonodomain=true)." << endl;
  }

  // Initialise external stimuli strictly from the monodomainStimulus
  // sub-dictionary.
  const MonodomainStimulusProtocol monodomainStimulus =
      stimulusIO::loadMonodomainStimulusProtocol(electroProperties());

  if (monodomainStimulus.boxes.size() > 0) {
    List<labelHashSet> stimCellSets(monodomainStimulus.boxes.size());
    labelHashSet stimCellSet;
    const vectorField &CI = mesh().C();
    forAll(CI, cellI) {
      const point &c = CI[cellI];
      forAll(monodomainStimulus.boxes, bI) {
        if (monodomainStimulus.boxes[bI].contains(c)) {
          stimCellSets[bI].insert(cellI);
          stimCellSet.insert(cellI);
        }
      }
    }

    // Set the number of stimulii
    externalStimulus_.setSize(monodomainStimulus.boxes.size());

    // Initialise each stimulii
    forAll(monodomainStimulus.boxes, bI) {
      externalStimulus_[bI].cellIDs_ = stimCellSets[bI].toc();
      externalStimulus_[bI].startTime_ = monodomainStimulus.startTimes[bI];
      externalStimulus_[bI].duration_ = monodomainStimulus.durations[bI];
      externalStimulus_[bI].intensity_ = monodomainStimulus.intensities[bI];

      if (externalStimulus_[bI].cellIDs_.size() == 0) {
        WarningInFunction << "Stimulus box index " << bI
                          << " does not contain any mesh cells." << nl;
      }
    }

    // Only print stimulus info if this is a 3D mesh
    if (mesh().nGeometricD() == 3) {
      Info << "---------------------------------------------" << nl
           << "Stimulus source dictionary: monodomainStimulus" << nl
           << "3D mesh detected — stimulus parameters:" << nl
           << "Stimulus boxes: " << monodomainStimulus.boxes.size() << nl
           << "Number of cells in stimulus region: " << stimCellSet.size() << nl
           << "Stimulus duration list: " << monodomainStimulus.durations << nl
           << "Stimulus start time(s): " << monodomainStimulus.startTimes << nl
           << "Stimulus intensity list: " << monodomainStimulus.intensities
           << nl << "---------------------------------------------" << endl;
    }
  } else {
    externalStimulus_.setSize(0);

    if (mesh().nGeometricD() == 3) {
      Info << "No monodomainStimulus sub-dictionary provided. "
           << "Defaulting to zero external stimulus current." << endl;
    }
  }

  // Write parameters
  Info << "Surface-to-volume ratio chi = " << chi_ << nl
       << "Membrane capacitance Cm = " << Cm_ << nl << endl;

  // Allocate model-declared output fields before lifecycle hooks run.
  const wordList names = ionicModelPtr_->exportedFieldNames();
  const wordList requiredPostProcessNames =
      prePostProcessor_->requiredPostProcessFieldNames(*ionicModelPtr_);
  outFields_.setSize(names.size());
  forAll(names, i) {
    outFields_.set(
        i, new volScalarField(IOobject(names[i], runTime.timeName(), mesh(),
                                       IOobject::NO_READ, IOobject::AUTO_WRITE),
                              mesh(), dimless, "zeroGradient"));
  }

  electroModelsPrePostProcessor::validatePostProcessFields(
      ionicModelPtr_->type(), names, requiredPostProcessNames);

  electroModelsPrePostProcessor::allocateFields(preProcessFieldNames_, mesh(),
                                                "preProcess_",
                                                preProcessFields_);

  if (hasCustomPrePostProcessor(prePostProcessor_())) {
    Info << "Using pre/post processor " << prePostProcessor_->type()
         << " for ionic model " << ionicModelPtr_->type() << "." << endl;
  }

  if (hasCustomPrePostProcessor(prePostProcessor_()) &&
      !preProcessFieldNames_.empty()) {
    Info << "Running preProcess with fields " << preProcessFieldNames_ << "."
         << endl;
  }

  prePostProcessor_->preProcess(*ionicModelPtr_, Vm_, preProcessFields_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void monoDomainElectro::setDeltaT(Time &runTime) {
  // Update the time step
  if (solutionAlg() == solutionAlgorithm::EXPLICIT && setDeltaT_) {
    setDeltaT_ = false;

    // Unit face normals
    surfaceVectorField n("n", mesh().Sf());
    n /= mesh().magSf();

    // Interpolate the conductivity to the faces and take the normal
    // component and divide by chi and Cm
    const scalarField Df((n & (n & fvc::interpolate(conductivity_))) /
                         (chi_ * Cm_));

    // Lookup the desired Courant number
    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

    // Use the face delta coeffs as measures of the local cell size
    const scalarField dx(1.0 / mesh().deltaCoeffs());

    // Compute stable dt for every face and take the minium
    scalar newDeltaT = maxCo * gMin(sqr(dx) / Df);

    Info << "Setting deltaT = " << newDeltaT << ", maxCo = " << maxCo << endl;

    runTime.setDeltaT(newDeltaT);
  }
}

bool monoDomainElectro::evolve() {
  bool converged = true;
  if (solutionAlg() == solutionAlgorithm::IMPLICIT) {
    converged = evolveImplicit();
  } else if (solutionAlg() == solutionAlgorithm::EXPLICIT) {
    converged = evolveExplicit();
  } else {
    FatalErrorInFunction
        << "Unrecognised solution algorithm. Available options are "
        << electroModel::solutionAlgorithmNames_
               [electroModel::solutionAlgorithm::IMPLICIT]
        << ", "
        << electroModel::solutionAlgorithmNames_
               [electroModel::solutionAlgorithm::EXPLICIT]
        << endl;
  }

  const bool shouldPostProcess =
      prePostProcessor_->shouldPostProcess(*ionicModelPtr_, Vm_);

  if ((runTime().outputTime() || shouldPostProcess) && !outFields_.empty()) {
    ionicModelPtr_->exportStates(outFields_);
  }

  if (shouldPostProcess) {
    if (hasCustomPrePostProcessor(prePostProcessor_())) {
      Info << "Running postProcess via " << prePostProcessor_->type() << "."
           << endl;
    }
    prePostProcessor_->postProcess(*ionicModelPtr_, Vm_, outFields_);
  }

  return converged;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
