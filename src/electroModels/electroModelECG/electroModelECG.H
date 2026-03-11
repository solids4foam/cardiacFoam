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

#ifndef electroModelECG_H
#define electroModelECG_H

#include "physicsModel.H"
#include "electroModel.H"
#include "ecgModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace electroModels {

/*---------------------------------------------------------------------------*\
                           Class electroModelECG Declaration
\*---------------------------------------------------------------------------*/

class electroModelECG : public physicsModel {
  //- The underlying electro model
  autoPtr<electroModel> electroModelPtr_;

  //- The ECG model
  autoPtr<ecgModel> ecgModelPtr_;

  //- Disallow default bitwise copy construct
  electroModelECG(const electroModelECG &);

  //- Disallow default bitwise assignment
  void operator=(const electroModelECG &);

public:
  //- Runtime type information
  TypeName("electroModel+ECG");

  //- Runtime selection constructor entrypoint
  static autoPtr<physicsModel> New(
      Time &runTime,
      const word &region = dynamicFvMesh::defaultRegion);

  // Constructors

  //- Construct from components
  electroModelECG(Time &runTime,
                  const word &region = dynamicFvMesh::defaultRegion);

  //- Construct from components
  electroModelECG(const word &type, Time &runTime,
                  const word &region = dynamicFvMesh::defaultRegion);

  //- Destructor
  virtual ~electroModelECG();

  // Member Functions

  // Edit

  //- Update the size of the time-step
  virtual void setDeltaT(Time &runTime);

  //- Evolve the electro model
  virtual bool evolve();

  //- Read electroProperties dictionary
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
