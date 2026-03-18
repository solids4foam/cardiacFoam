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

Class
    Foam::manufacturedFDA

Description
    Wrapper class for calling from the cellML ionic model for the ten Tusscher,
    Noble, Noble, Panfilov, 2004, model.

SourceFiles
    manufacturedFDA.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef manufacturedFDA_H
#define manufacturedFDA_H

#include "ionicModel.H"
#include "tmanufacturedFDA_2014Names.H"
#include "OFstream.H"
#include "volFields.H"

extern const char *tmanufacturedFDASTATES_NAMES[NUM_STATES];
extern const char *tmanufacturedFDAALGEBRAIC_NAMES[NUM_ALGEBRAIC];
extern const char *tmanufacturedFDACONSTANTS_NAMES[NUM_CONSTANTS];

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                         Class manufacturedFDA Declaration
\*---------------------------------------------------------------------------*/

class tmanufacturedFDA : public ionicModel {
  // Private Data

  //- States
  mutable PtrList<scalarField> STATES_;

  //- Constants
  mutable scalarField CONSTANTS_;

  //- Algebraic
  mutable PtrList<scalarField> ALGEBRAIC_;

  //- Rates
  mutable PtrList<scalarField> RATES_;

  // Private Member Functions

  //- No copy construct
  tmanufacturedFDA(const tmanufacturedFDA &) = delete;

  //- No copy assignment
  void operator=(const tmanufacturedFDA &) = delete;

protected:
  virtual const PtrList<scalarField> *ioStatesPtr() const override {
    return &STATES_;
  }

  virtual const PtrList<scalarField> *ioAlgebraicPtr() const override {
    return &ALGEBRAIC_;
  }

  virtual const PtrList<scalarField> *ioRatesPtr() const override {
    return &RATES_;
  }

  virtual const scalarField *ioConstantsPtr() const override {
    return &CONSTANTS_;
  }

  virtual label ioNumConstants() const override { return NUM_CONSTANTS; }

  virtual const char *const *ioConstantNames() const override {
    return tmanufacturedFDACONSTANTS_NAMES;
  }

  virtual const char *const *ioStateNames() const override {
    return tmanufacturedFDASTATES_NAMES;
  }

  virtual label ioNumStates() const override { return NUM_STATES; }

  virtual const char *const *ioAlgebraicNames() const override {
    return tmanufacturedFDAALGEBRAIC_NAMES;
  }

  virtual label ioNumAlgebraic() const override { return NUM_ALGEBRAIC; }

public:
  //- Runtime type information
  TypeName("tmanufacturedFDA");

  // Override the only available dimension family in this model
  virtual List<word> supportedDimensions() const override;

  virtual word tissueName() const override
  {
    return (tissue() == 1) ? "1D" :
           (tissue() == 2) ? "2D" :
           (tissue() == 3) ? "3D" :
                             "unknown";
  }

  // Constructors

  //- Construct given a dictionary, the number of integration points and a
  //  initial time step for the ODE solver
  tmanufacturedFDA(const dictionary &dict, const label nIntegrationPoints,
                   const scalar initialDeltaT,
                   const Switch solveVmWithinODESolver = false);

  //- Destructor
  virtual ~tmanufacturedFDA();

  // Member Functions

  virtual void solveODE(const scalar stepStartTime, const scalar deltaT,
                        const scalarField &Vm, scalarField &Im) override;

  //- Number of ODE's to solve
  virtual label nEqns() const override { return NUM_STATES; }

  //- Derivatives of each ODE with respect to the primary unknown in that
  //  equation
  virtual void derivatives(const scalar t, const scalarField &y,
                           scalarField &dydt) const override;

  // to provide variable names for export
  virtual Foam::wordList exportedFieldNames() const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
