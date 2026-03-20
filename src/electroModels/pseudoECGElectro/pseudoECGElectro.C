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
#include "DynamicList.H"
#include "ecgModelIO.H"
#include "pseudoECGManufacturedReference.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pseudoECGElectro, 0);
addToRunTimeSelectionTable(ecgModel, pseudoECGElectro, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

pseudoECGElectro::pseudoECGElectro
(
    const volScalarField& Vm,
    const dictionary& dict,
    const volTensorField& conductivity
)
: ecgModel(Vm, dict, conductivity), outputPtr_(), manufacturedOutputPtr_(),
      manufacturedEnabled_(false),
      manufacturedDimension_(max(label(1), min(mesh_.nGeometricD(), label(3)))),
      manufacturedReferenceQuadratureOrder_(96),
      manufacturedCheckQuadratureOrders_(),
      manufacturedReferenceSpatialValues_(electrodePositions_.size(), scalar(0)),
      manufacturedCheckSpatialValues_(),
      manufacturedCachedConductivity_(tensor::zero),
      manufacturedSpatialCacheValid_(false), manufacturedSampleCount_(0),
      referenceErrorL1Sum_(electrodePositions_.size(), scalar(0)),
      referenceErrorL2Sum_(electrodePositions_.size(), scalar(0)),
      referenceErrorLinf_(electrodePositions_.size(), scalar(0)),
      checkErrorL1Sum_(),
      checkErrorL2Sum_(),
      checkErrorLinf_(),
      referenceDeltaL1Sum_(),
      referenceDeltaL2Sum_(),
      referenceDeltaLinf_(),
      manufacturedSummaryWritten_(false)
{
  const fileName outDir(mesh_.time().path() / "postProcessing");
  outputPtr_ =
      ecgModelIO::openTimeSeries(outDir, "pseudoECG.dat", electrodeNames_);
  readManufacturedConfig(dict);
  if (manufacturedEnabled_)
  {
    initialiseManufacturedOutput();
  }
}

void pseudoECGElectro::resizeManufacturedCheckStorage()
{
  const label nChecks = manufacturedCheckQuadratureOrders_.size();
  const label nElectrodes = electrodePositions_.size();

  manufacturedCheckSpatialValues_.setSize(nChecks);
  checkErrorL1Sum_.setSize(nChecks);
  checkErrorL2Sum_.setSize(nChecks);
  checkErrorLinf_.setSize(nChecks);
  referenceDeltaL1Sum_.setSize(nChecks);
  referenceDeltaL2Sum_.setSize(nChecks);
  referenceDeltaLinf_.setSize(nChecks);

  for (label checkI = 0; checkI < nChecks; ++checkI) {
    manufacturedCheckSpatialValues_[checkI].setSize(nElectrodes, scalar(0));
    checkErrorL1Sum_[checkI].setSize(nElectrodes, scalar(0));
    checkErrorL2Sum_[checkI].setSize(nElectrodes, scalar(0));
    checkErrorLinf_[checkI].setSize(nElectrodes, scalar(0));
    referenceDeltaL1Sum_[checkI].setSize(nElectrodes, scalar(0));
    referenceDeltaL2Sum_[checkI].setSize(nElectrodes, scalar(0));
    referenceDeltaLinf_[checkI].setSize(nElectrodes, scalar(0));
  }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pseudoECGElectro::readManufacturedConfig(const dictionary& dict)
{
  manufacturedEnabled_ = false;
  manufacturedDimension_ = max(label(1), min(mesh_.nGeometricD(), label(3)));
  manufacturedReferenceQuadratureOrder_ = 96;
  manufacturedCheckQuadratureOrders_.setSize(4);
  manufacturedCheckQuadratureOrders_[0] = 6;
  manufacturedCheckQuadratureOrders_[1] = 12;
  manufacturedCheckQuadratureOrders_[2] = 24;
  manufacturedCheckQuadratureOrders_[3] = 48;

  if (!dict.found("manufactured"))
  {
    resizeManufacturedCheckStorage();
    return;
  }

  const dictionary &manufacturedDict = dict.subDict("manufactured");
  manufacturedEnabled_ =
      manufacturedDict.lookupOrDefault<Switch>("enabled", true);

  if (!manufacturedEnabled_)
  {
    resizeManufacturedCheckStorage();
    return;
  }

  if (manufacturedDict.found("dimension"))
  {
    const word dimensionName(manufacturedDict.lookup("dimension"));
    if (dimensionName == "1D")
    {
      manufacturedDimension_ = 1;
    }
    else if (dimensionName == "2D")
    {
      manufacturedDimension_ = 2;
    }
    else if (dimensionName == "3D")
    {
      manufacturedDimension_ = 3;
    }
    else
    {
      FatalErrorInFunction
          << "Unsupported manufactured pseudo-ECG dimension '"
          << dimensionName << "'. Expected one of 1D, 2D, or 3D."
          << exit(FatalError);
    }
  }

  manufacturedReferenceQuadratureOrder_ =
      manufacturedDict.lookupOrDefault<label>("referenceQuadratureOrder", 96);
  if (manufacturedDict.found("checkQuadratureOrders"))
  {
    manufacturedDict.lookup("checkQuadratureOrders")
        >> manufacturedCheckQuadratureOrders_;
  }
  else
  {
    manufacturedCheckQuadratureOrders_.setSize(1);
    manufacturedCheckQuadratureOrders_[0] =
        manufacturedDict.lookupOrDefault<label>("checkQuadratureOrder", 6);
  }

  if (manufacturedCheckQuadratureOrders_.size() == 0)
  {
    FatalErrorInFunction
        << "Manufactured pseudo-ECG requires at least one entry in "
           "checkQuadratureOrders."
        << exit(FatalError);
  }

  for (label i = 0; i < manufacturedCheckQuadratureOrders_.size(); ++i)
  {
    for (label j = i + 1; j < manufacturedCheckQuadratureOrders_.size(); ++j)
    {
      if (manufacturedCheckQuadratureOrders_[j] <
          manufacturedCheckQuadratureOrders_[i]) {
        const label tmp = manufacturedCheckQuadratureOrders_[i];
        manufacturedCheckQuadratureOrders_[i] =
            manufacturedCheckQuadratureOrders_[j];
        manufacturedCheckQuadratureOrders_[j] = tmp;
      }
    }
  }
  label uniqueCount = 0;
  forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
    const label current = manufacturedCheckQuadratureOrders_[checkI];
    if (checkI == 0 ||
        current != manufacturedCheckQuadratureOrders_[uniqueCount - 1]) {
      manufacturedCheckQuadratureOrders_[uniqueCount++] = current;
    }
  }
  manufacturedCheckQuadratureOrders_.setSize(uniqueCount);

  if (!pseudoECGManufacturedSupportsQuadratureOrder(
          manufacturedReferenceQuadratureOrder_)) {
    FatalErrorInFunction
        << "Manufactured pseudo-ECG quadrature orders must be positive. "
        << "Requested qRef=" << manufacturedReferenceQuadratureOrder_ << "."
        << exit(FatalError);
  }

  forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
    if (!pseudoECGManufacturedSupportsQuadratureOrder(
            manufacturedCheckQuadratureOrders_[checkI])) {
      FatalErrorInFunction
          << "Manufactured pseudo-ECG quadrature orders must be positive. "
          << "Invalid qCheck=" << manufacturedCheckQuadratureOrders_[checkI]
          << "."
          << exit(FatalError);
    }
  }

  resizeManufacturedCheckStorage();

  Info << "Enabled manufactured pseudo-ECG references: qRef="
       << manufacturedReferenceQuadratureOrder_ << ", qChecks=(";
  forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
    if (checkI > 0)
    {
      Info << " ";
    }
    Info << manufacturedCheckQuadratureOrders_[checkI];
  }
  Info << "), dimension=" << manufacturedDimension_ << "." << endl;
}

void pseudoECGElectro::initialiseManufacturedOutput()
{
  const fileName outDir(mesh_.time().path() / "postProcessing");
  wordList columns;

  forAll(electrodeNames_, electrodeI)
  {
    const word &name = electrodeNames_[electrodeI];
    columns.append("numeric_" + name);
    forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
      columns.append(
          "refQ" + Foam::name(manufacturedCheckQuadratureOrders_[checkI]) + "_" +
          name);
    }
    columns.append("refQ" +
                   Foam::name(manufacturedReferenceQuadratureOrder_) + "_" +
                   name);
    forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
      columns.append(
          "errQ" + Foam::name(manufacturedCheckQuadratureOrders_[checkI]) + "_" +
          name);
    }
    columns.append("errQ" +
                   Foam::name(manufacturedReferenceQuadratureOrder_) + "_" +
                   name);
    forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
      columns.append(
          "deltaQuadratureQ" +
          Foam::name(manufacturedCheckQuadratureOrders_[checkI]) + "_Q" +
          Foam::name(manufacturedReferenceQuadratureOrder_) + "_" + name);
    }
  }

  manufacturedOutputPtr_ =
      ecgModelIO::openTimeSeries(outDir, "manufacturedPseudoECG.dat", columns);
}

void pseudoECGElectro::invalidateManufacturedReferenceCache()
{
  manufacturedReferenceSpatialValues_.clear();
  manufacturedCheckSpatialValues_.clear();
  manufacturedSpatialCacheValid_ = false;
  manufacturedCachedConductivity_ = tensor::zero;
}

void pseudoECGElectro::rebuildManufacturedReferenceCache
(
    const tensor& referenceConductivity
)
{
  List<scalar> referenceNodes, referenceWeights;
  pseudoECGManufacturedQuadratureRule(manufacturedReferenceQuadratureOrder_,
                                      referenceNodes, referenceWeights);
  List<List<scalar>> checkNodes(manufacturedCheckQuadratureOrders_.size());
  List<List<scalar>> checkWeights(manufacturedCheckQuadratureOrders_.size());
  forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
    pseudoECGManufacturedQuadratureRule(
        manufacturedCheckQuadratureOrders_[checkI], checkNodes[checkI],
        checkWeights[checkI]);
  }

  manufacturedReferenceSpatialValues_.setSize(electrodePositions_.size());
  manufacturedCheckSpatialValues_.setSize(manufacturedCheckQuadratureOrders_.size());
  forAll(manufacturedCheckSpatialValues_, checkI)
  {
    manufacturedCheckSpatialValues_[checkI].setSize(electrodePositions_.size());
  }

  forAll(electrodePositions_, electrodeI)
  {
    forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
      manufacturedCheckSpatialValues_[checkI][electrodeI] =
          computeManufacturedPseudoECGSpatialReference(
              electrodePositions_[electrodeI], referenceConductivity,
              manufacturedDimension_, checkNodes[checkI], checkWeights[checkI]);
    }
    manufacturedReferenceSpatialValues_[electrodeI] =
        computeManufacturedPseudoECGSpatialReference(
            electrodePositions_[electrodeI], referenceConductivity,
            manufacturedDimension_, referenceNodes, referenceWeights);
  }

  manufacturedCachedConductivity_ = referenceConductivity;
  manufacturedSpatialCacheValid_ = true;
}

void pseudoECGElectro::computePseudoECGValues(List<scalar>& values) const
{
  // Gima-Rudy dipole:
  //   phi_pseudo(P) = -sum_c [ (conductivity . grad(Vm))_c . r_vec * V_c / |r|^3 ]
  //   r_vec = C_c - P

  const tmp<volVectorField> tgradVm = fvc::grad(Vm_);
  const vectorField &gradVm = tgradVm().primitiveField();

  const scalarField &Vols = mesh_.V();
  const vectorField &Ctrs = mesh_.C().primitiveField();
  const tensorField &conductivityField = conductivity().primitiveField();

  const label nE = electrodePositions_.size();

  values.setSize(nE);
  forAll(values, valueI)
  {
    values[valueI] = scalar(0);
  }

  forAll(Ctrs, cI)
  {
    const vector dipole = (conductivityField[cI] & gradVm[cI]) * Vols[cI];

    for (label pI = 0; pI < nE; pI++)
    {
      const vector r_vec = Ctrs[cI] - electrodePositions_[pI];
      const scalar r = mag(r_vec);
      if (r > VSMALL)
      {
        values[pI] += (dipole & r_vec) / (r * r * r);
      }
    }
  }

  // Parallel reduction (no sign flip: consistent with original accumulation)
  for (label pI = 0; pI < nE; pI++)
  {
    reduce(values[pI], sumOp<scalar>());
  }
}

void pseudoECGElectro::updateManufacturedStatistics
(
    const List<scalar>& numericValues,
    const List<List<scalar>>& checkReferenceValues,
    const List<scalar>& referenceValues
)
{
  ++manufacturedSampleCount_;

  forAll(numericValues, electrodeI) {
    const scalar referenceError =
        Foam::mag(numericValues[electrodeI] - referenceValues[electrodeI]);

    referenceErrorL1Sum_[electrodeI] += referenceError;
    referenceErrorL2Sum_[electrodeI] += referenceError * referenceError;
    referenceErrorLinf_[electrodeI] =
        max(referenceErrorLinf_[electrodeI], referenceError);
    forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
      const scalar checkError = Foam::mag(
          numericValues[electrodeI] - checkReferenceValues[checkI][electrodeI]);
      const scalar referenceDelta =
          Foam::mag(referenceValues[electrodeI] -
                    checkReferenceValues[checkI][electrodeI]);

      checkErrorL1Sum_[checkI][electrodeI] += checkError;
      checkErrorL2Sum_[checkI][electrodeI] += checkError * checkError;
      checkErrorLinf_[checkI][electrodeI] =
          max(checkErrorLinf_[checkI][electrodeI], checkError);

      referenceDeltaL1Sum_[checkI][electrodeI] += referenceDelta;
      referenceDeltaL2Sum_[checkI][electrodeI] += referenceDelta * referenceDelta;
      referenceDeltaLinf_[checkI][electrodeI] =
          max(referenceDeltaLinf_[checkI][electrodeI], referenceDelta);
    }
  }
}

void pseudoECGElectro::writeManufacturedSummary()
{
  if (!manufacturedEnabled_ || manufacturedSummaryWritten_)
  {
    return;
  }

  manufacturedSummaryWritten_ = true;

  if (!Pstream::master())
  {
    return;
  }

  const fileName outputFile(mesh_.time().path() / "postProcessing" /
                            "manufacturedPseudoECGSummary.dat");
  OFstream os(outputFile);

  os << "Manufactured pseudo-ECG summary\n";
  os << "samples " << manufacturedSampleCount_ << "\n";
  os << "dimension " << manufacturedDimension_ << "D\n";
  os << "qChecks";
  forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
    os << " " << manufacturedCheckQuadratureOrders_[checkI];
  }
  os << "\n";
  os << "qReference " << manufacturedReferenceQuadratureOrder_ << "\n";
  os << "Electrode  L1_err_ref  L2_err_ref  Linf_err_ref";
  forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
    const label qCheck = manufacturedCheckQuadratureOrders_[checkI];
    os << "  L1_err_q" << qCheck << "  L2_err_q" << qCheck << "  Linf_err_q"
       << qCheck << "  L1_delta_q" << qCheck << "_ref  L2_delta_q" << qCheck
       << "_ref  Linf_delta_q" << qCheck << "_ref";
  }
  os << "\n";

  const scalar sampleCount = max(scalar(1), scalar(manufacturedSampleCount_));

  forAll(electrodeNames_, electrodeI)
  {
    os << electrodeNames_[electrodeI] << " "
       << referenceErrorL1Sum_[electrodeI] / sampleCount << " "
       << Foam::sqrt(referenceErrorL2Sum_[electrodeI] / sampleCount) << " "
       << referenceErrorLinf_[electrodeI];
    forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
      os << " " << checkErrorL1Sum_[checkI][electrodeI] / sampleCount << " "
         << Foam::sqrt(checkErrorL2Sum_[checkI][electrodeI] / sampleCount)
         << " " << checkErrorLinf_[checkI][electrodeI] << " "
         << referenceDeltaL1Sum_[checkI][electrodeI] / sampleCount << " "
         << Foam::sqrt(referenceDeltaL2Sum_[checkI][electrodeI] / sampleCount)
         << " " << referenceDeltaLinf_[checkI][electrodeI];
    }
    os << "\n";
  }
}

void pseudoECGElectro::compute()
{
  List<scalar> numericValues;
  computePseudoECGValues(numericValues);

  if (manufacturedEnabled_)
  {
    const scalar timeValue = mesh_.time().value();
    const tensor referenceConductivity = conductivity().primitiveField()[0];
    if (!manufacturedSpatialCacheValid_ ||
        Foam::magSqr(referenceConductivity - manufacturedCachedConductivity_) >
            SMALL) {
      rebuildManufacturedReferenceCache(referenceConductivity);
    }

    const scalar timeFactor = pseudoECGManufacturedTimeFactor(timeValue);
    List<List<scalar>> checkReferenceValues(
        manufacturedCheckQuadratureOrders_.size());
    forAll(checkReferenceValues, checkI) {
      checkReferenceValues[checkI].setSize(electrodePositions_.size(),
                                           scalar(0));
    }
    List<scalar> referenceValues(electrodePositions_.size(), scalar(0));
    DynamicList<scalar> manufacturedRow;
    manufacturedRow.reserve(
        electrodePositions_.size() *
        (2 + 3 * manufacturedCheckQuadratureOrders_.size()));

    forAll(electrodePositions_, electrodeI)
  {
      referenceValues[electrodeI] =
          timeFactor * manufacturedReferenceSpatialValues_[electrodeI];
      manufacturedRow.append(numericValues[electrodeI]);
      forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
        checkReferenceValues[checkI][electrodeI] =
            timeFactor * manufacturedCheckSpatialValues_[checkI][electrodeI];
        manufacturedRow.append(checkReferenceValues[checkI][electrodeI]);
      }
      manufacturedRow.append(referenceValues[electrodeI]);
      forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
        manufacturedRow.append(Foam::mag(
            numericValues[electrodeI] - checkReferenceValues[checkI][electrodeI]));
      }
      manufacturedRow.append(
          Foam::mag(numericValues[electrodeI] - referenceValues[electrodeI]));
      forAll(manufacturedCheckQuadratureOrders_, checkI)
  {
        manufacturedRow.append(
            Foam::mag(referenceValues[electrodeI] -
                      checkReferenceValues[checkI][electrodeI]));
      }
    }

    updateManufacturedStatistics(numericValues, checkReferenceValues,
                                 referenceValues);
    List<scalar> manufacturedValues(manufacturedRow.size(), scalar(0));
    forAll(manufacturedValues, valueI) {
      manufacturedValues[valueI] = manufacturedRow[valueI];
    }
    ecgModelIO::writeRow(manufacturedOutputPtr_.ref(), timeValue,
                         manufacturedValues);

    const Time &time = mesh_.time();
    if (time.value() + 0.5 * time.deltaTValue() >= time.endTime().value()) {
      writeManufacturedSummary();
    }
  }

  if (mesh_.time().outputTime())
  {
    ecgModelIO::writeRow(outputPtr_.ref(), mesh_.time().value(), numericValues);
  }
}

bool pseudoECGElectro::read(const dictionary& dict)
{
  const wordList previousElectrodeNames(electrodeNames_);
  const List<vector> previousElectrodePositions(electrodePositions_);
  const bool wasManufacturedEnabled = manufacturedEnabled_;
  const label previousManufacturedDimension = manufacturedDimension_;
  const label previousReferenceQuadratureOrder =
      manufacturedReferenceQuadratureOrder_;
  const List<label> previousCheckQuadratureOrders =
      manufacturedCheckQuadratureOrders_;

  if (!ecgModel::read(dict)) {
    return false;
  }

  if (electrodeNames_.size() != previousElectrodeNames.size()) {
    FatalErrorInFunction
        << "Changing the number of ECG electrodes during read() is not "
           "supported because the output column layout is fixed at startup."
        << exit(FatalError);
  }

  forAll(previousElectrodeNames, electrodeI) {
    if (electrodeNames_[electrodeI] != previousElectrodeNames[electrodeI]) {
      FatalErrorInFunction
          << "Changing ECG electrode names during read() is not supported "
             "because the output column layout is fixed at startup."
          << exit(FatalError);
    }
  }

  bool electrodePositionsChanged = false;
  forAll(previousElectrodePositions, electrodeI) {
    if (Foam::magSqr(electrodePositions_[electrodeI] -
                     previousElectrodePositions[electrodeI]) > SMALL) {
      electrodePositionsChanged = true;
      break;
    }
  }

  readManufacturedConfig(dict);

  const bool manufacturedConfigurationChanged =
      manufacturedEnabled_ != wasManufacturedEnabled ||
      manufacturedDimension_ != previousManufacturedDimension ||
      manufacturedReferenceQuadratureOrder_ != previousReferenceQuadratureOrder ||
      manufacturedCheckQuadratureOrders_ != previousCheckQuadratureOrders;

  if (manufacturedOutputPtr_.valid() && manufacturedConfigurationChanged &&
      manufacturedEnabled_) {
    FatalErrorInFunction
        << "Changing manufactured pseudo-ECG dimension or quadrature orders "
           "during read() is not supported because the manufactured output "
           "layout and semantics are fixed at startup."
        << exit(FatalError);
  }

  if (!manufacturedEnabled_)
  {
    invalidateManufacturedReferenceCache();
  } else if (electrodePositionsChanged || manufacturedConfigurationChanged) {
    invalidateManufacturedReferenceCache();
  }

  if (manufacturedEnabled_ && !wasManufacturedEnabled &&
      !manufacturedOutputPtr_.valid()) {
    initialiseManufacturedOutput();
  }

  return true;
}

void pseudoECGElectro::end()
{
  writeManufacturedSummary();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
