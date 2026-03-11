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

#include "electroModelECG.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace electroModels {
defineTypeNameAndDebug(electroModelECG, 0);
addToRunTimeSelectionTable(physicsModel, electroModelECG, physicsModel);
} // namespace electroModels
} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electroModels::electroModelECG::electroModelECG(const word &type,
                                                      Time &runTime,
                                                      const word &region)
    : physicsModel(type, runTime),
      electroModelPtr_(electroModel::New(runTime, region)), ecgModelPtr_() {
  const dictionary &props = electroModelPtr_->electroProperties();

  if (props.found("ECG")) {
    const dictionary &ecgDict = props.subDict("ECG");
    ecgModelPtr_ = ecgModel::New(electroModelPtr_->Vm(), ecgDict);
  } else {
    FatalErrorInFunction
        << "ECG sub-dictionary not found in electroProperties for "
        << "electroModel+ECG type." << abort(FatalError);
  }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroModels::electroModelECG::~electroModelECG() {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electroModels::electroModelECG::setDeltaT(Time &runTime) {
  electroModelPtr_->setDeltaT(runTime);
}

bool Foam::electroModels::electroModelECG::evolve() {
  bool status = electroModelPtr_->evolve();

  if (ecgModelPtr_.valid()) {
    ecgModelPtr_->compute();
  }

  return status;
}

bool Foam::electroModels::electroModelECG::read() {
  bool status = electroModelPtr_->read();

  if (ecgModelPtr_.valid()) {
    const dictionary &props = electroModelPtr_->electroProperties();
    if (props.found("ECG")) {
      ecgModelPtr_->read(props.subDict("ECG"));
    }
  }

  return status;
}

void Foam::electroModels::electroModelECG::writeFields(const Time &runTime) {
  electroModelPtr_->writeFields(runTime);
}

void Foam::electroModels::electroModelECG::end() { electroModelPtr_->end(); }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
