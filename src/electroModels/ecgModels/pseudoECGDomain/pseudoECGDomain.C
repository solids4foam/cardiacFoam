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

#include "pseudoECGDomain.H"
#include "addToRunTimeSelectionTable.H"
#include "ecgModelIO.H"
#include "ecgVerificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PseudoECGDomain, 0);
addToRunTimeSelectionTable(ECGDomain, PseudoECGDomain, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace
{

void finalizeVerificationModel(autoPtr<ecgVerificationModel>& verifierPtr)
{
    if (verifierPtr.valid())
    {
        verifierPtr->end();
        verifierPtr.clear();
    }
}

}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

PseudoECGDomain::PseudoECGDomain
(
    const electroStateProvider& stateProvider,
    const dictionary& dict,
    const word& domainName
)
: ECGDomain(stateProvider, dict, domainName),
  outputPtr_(),
  solverPtr_(ECGSolver::New(dict)),
  verificationModelPtr_()
{
  const fileName outDir(mesh_.time().path() / "postProcessing");
  outputPtr_ =
      ecgModelIO::openTimeSeries(outDir, "pseudoECG.dat", electrodeNames_);
  verificationModelPtr_ = ecgVerificationModel::New
  (
      stateProvider,
      dict,
      electrodeNames_,
      electrodePositions_
  );
  numericValues_.setSize(electrodeNames_.size(), 0.0);
}

PseudoECGDomain::~PseudoECGDomain() = default;

void PseudoECGDomain::evolve
(
    scalar t0,
    scalar dt
)
{
  (void)t0;
  (void)dt;

  solverPtr_->compute(*this, numericValues_);

  if (verificationModelPtr_.valid())
  {
    verificationModelPtr_->record(numericValues_);
  }

  if (mesh_.time().outputTime())
  {
    ecgModelIO::writeRow
    (
        outputPtr_.ref(), mesh_.time().value(), numericValues_
    );
  }
}

bool PseudoECGDomain::read(const dictionary& dict)
{
  const wordList previousElectrodeNames(electrodeNames_);

  if (!ECGDomain::read(dict)) {
    return false;
  }

  numericValues_.setSize(electrodeNames_.size(), 0.0);

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

  const word requestedType(ecgVerificationModel::selectedType(dict));

  if (requestedType.empty())
  {
    finalizeVerificationModel(verificationModelPtr_);
    return true;
  }

  if
  (
      verificationModelPtr_.valid()
   && verificationModelPtr_->type() == requestedType
  )
  {
    verificationModelPtr_->updateElectrodes
    (
        electrodeNames_,
        electrodePositions_
    );

    return verificationModelPtr_->read(dict);
  }

  finalizeVerificationModel(verificationModelPtr_);
  verificationModelPtr_ = ecgVerificationModel::New
  (
      stateProvider_,
      dict,
      electrodeNames_,
      electrodePositions_
  );

  return true;
}

void PseudoECGDomain::end()
{
  if (verificationModelPtr_.valid())
  {
    verificationModelPtr_->end();
  }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
