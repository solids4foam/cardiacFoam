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

#ifndef electroModelsPrePostProcessor_H
#define electroModelsPrePostProcessor_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "volFields.H"

namespace Foam {

class ionicModel;

namespace electroModels {

class electroModelsPrePostProcessor {
  dictionary dict_;

public:
  TypeName("electroModelsPrePostProcessor");

  declareRunTimeSelectionTable(autoPtr, electroModelsPrePostProcessor, dictionary,
                               (const dictionary &dict), (dict));

  static autoPtr<electroModelsPrePostProcessor> New(const word &ionicModelType,
                                                    const dictionary &dict);

  explicit electroModelsPrePostProcessor(const dictionary &dict);

  virtual ~electroModelsPrePostProcessor();

  static void allocateFields(const wordList &names, const fvMesh &mesh,
                             const word &prefix,
                             PtrList<volScalarField> &fields);

  static void validatePostProcessFields(const word &ionicModelType,
                                        const wordList &availableNames,
                                        const wordList &requiredNames);

  virtual void preProcess(ionicModel &model, volScalarField &Vm,
                          PtrList<volScalarField> &outFields);

  virtual wordList preProcessFieldNames(const ionicModel &model) const;

  virtual wordList requiredPostProcessFieldNames(const ionicModel &model) const;

  virtual bool shouldPostProcess(const ionicModel &model,
                                 const volScalarField &Vm) const;

  virtual void postProcess(const ionicModel &model, const volScalarField &Vm,
                           const PtrList<volScalarField> &fields);

protected:
  const dictionary &dict() const { return dict_; }
};

} // End namespace electroModels
} // End namespace Foam

#endif
