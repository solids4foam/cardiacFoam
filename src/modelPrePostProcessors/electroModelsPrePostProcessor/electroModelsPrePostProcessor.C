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

#include "electroModelsPrePostProcessor.H"

namespace Foam {
namespace electroModels {

defineTypeNameAndDebug(electroModelsPrePostProcessor, 0);
defineRunTimeSelectionTable(electroModelsPrePostProcessor, dictionary);

electroModelsPrePostProcessor::electroModelsPrePostProcessor(const dictionary &dict)
    : dict_(dict) {}

electroModelsPrePostProcessor::~electroModelsPrePostProcessor() {}

void electroModelsPrePostProcessor::allocateFields(
    const wordList &names, const fvMesh &mesh, const word &prefix,
    PtrList<volScalarField> &fields) {
  fields.setSize(names.size());

  forAll(names, i) {
    fields.set(
        i, new volScalarField(IOobject(prefix + names[i], mesh.time().timeName(),
                                       mesh, IOobject::NO_READ,
                                       IOobject::NO_WRITE),
                              mesh, dimless, "zeroGradient"));
  }
}

void electroModelsPrePostProcessor::validatePostProcessFields(
    const word &ionicModelType, const wordList &availableNames,
    const wordList &requiredNames) {
  forAll(requiredNames, i) {
    if (!availableNames.found(requiredNames[i])) {
      FatalErrorInFunction
          << "Pre/post processor for ionic model " << ionicModelType
          << " requires exported field '" << requiredNames[i]
          << "' for post-processing, but the configured ionic export list is "
          << availableNames << "." << nl
          << "Add it under outputVariables/ionic/export."
          << exit(FatalError);
    }
  }
}

autoPtr<electroModelsPrePostProcessor> electroModelsPrePostProcessor::New(
    const word &ionicModelType, const dictionary &dict) {
  auto *ctorPtr = dictionaryConstructorTable(ionicModelType);

  if (ctorPtr) {
    return autoPtr<electroModelsPrePostProcessor>(ctorPtr(dict));
  }

  return autoPtr<electroModelsPrePostProcessor>(
      new electroModelsPrePostProcessor(dict));
}

void electroModelsPrePostProcessor::preProcess(ionicModel &, volScalarField &,
                                               PtrList<volScalarField> &) {}

wordList electroModelsPrePostProcessor::preProcessFieldNames(
    const ionicModel &) const {
  return wordList();
}

wordList
electroModelsPrePostProcessor::requiredPostProcessFieldNames(
    const ionicModel &) const {
  return wordList();
}

bool electroModelsPrePostProcessor::shouldPostProcess(
    const ionicModel &, const volScalarField &) const {
  return false;
}

void electroModelsPrePostProcessor::postProcess(
    const ionicModel &, const volScalarField &,
    const PtrList<volScalarField> &) {}

} // End namespace electroModels
} // End namespace Foam
