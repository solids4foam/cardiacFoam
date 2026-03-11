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

#ifndef electroMechanicalModel_H
#define electroMechanicalModel_H

#include "physicsModel.H"
#include "electroModel.H"
#include "activeTensionModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace electroModels {

/*---------------------------------------------------------------------------*\
                         Class electroMechanicalModel Declaration
\*---------------------------------------------------------------------------*/

class electroMechanicalModel : public physicsModel {
  //- The underlying electro model
  autoPtr<electroModel> electroModelPtr_;

  //- The Active Tension model
  autoPtr<activeTensionModel> activeTensionModelPtr_;

  //- Active tension field
  autoPtr<volScalarField> TaPtr_;

  //- Lambda field (stretch) - typically read from solid simulation
  //  For now we use a local one or access it if available.
  autoPtr<volScalarField> lambdaPtr_;

  //- Exported active-tension fields (states/algebraic/rates)
  PtrList<volScalarField> activeTensionOutFields_;

  //- Single-cell active-tension trace output (if mesh has one cell)
  autoPtr<OFstream> activeTensionOutputPtr_;

  //- Disallow default bitwise copy construct
  electroMechanicalModel(const electroMechanicalModel &);

  //- Disallow default bitwise assignment
  void operator=(const electroMechanicalModel &);

public:
  //- Runtime type information
  TypeName("electroMechanicalModel");

  //- Runtime selection constructor entrypoint
  static autoPtr<physicsModel> New(
      Time &runTime,
      const word &region = dynamicFvMesh::defaultRegion);

  // Constructors

  //- Construct from components
  electroMechanicalModel(
      Time &runTime,
      const word &region = dynamicFvMesh::defaultRegion);

  //- Construct from components
  electroMechanicalModel(const word &type, Time &runTime,
                         const word &region = dynamicFvMesh::defaultRegion);

  //- Destructor
  virtual ~electroMechanicalModel();

  // Member Functions

  // Edit

  //- Update the size of the time-step
  virtual void setDeltaT(Time &runTime);

  //- Evolve the models
  virtual bool evolve();

  //- Read properties
  virtual bool read();

  //- Write fields
  virtual void writeFields(const Time &runTime);

  //- This function is called at the end of the simulation
  virtual void end();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
