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
    Foam::pseudoECGElectro

Description
    ECG model using the Gima-Rudy dipole approximation (pseudo-ECG).

    The potential at each electrode position P is:

        phi_pseudo(P,t) = -sum_c [ (Gi . grad(Vm))_c . r_vec * V_c / |r|^3 ]

    where r_vec = C_c - P and the sum is made parallel-safe via Foam::reduce.
    No sigmaT factor is applied (consistent with original implementation).

    Dictionary (inside pseudoECGElectroCoeffs):

        pseudoECGElectroCoeffs
        {
            Gi      [-1 -3 3 0 0 2 0] (0.17 0 0  0.019 0  0.019); // fallback
            sigmaT  [-1 -3 3 0 0 2 0] 0.24725;   // stored but unused

            electrodes
            {
                V1  (-26.1184 -283.543 -70.0558);
                // ...
            }
        }

    Output: postProcessing/pseudoECG.dat  (at outputTime, one row per step)

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    pseudoECGElectro.C
\*---------------------------------------------------------------------------*/

#ifndef pseudoECGElectro_H
#define pseudoECGElectro_H

#include "ecgModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                     Class pseudoECGElectro Declaration
\*---------------------------------------------------------------------------*/

class pseudoECGElectro : public ecgModel {
  // Private data

  //- Output stream for pseudoECG.dat
  autoPtr<OFstream> outputPtr_;

  // Private Member Functions

  //- Disallow copy construct
  pseudoECGElectro(const pseudoECGElectro &);

  //- Disallow assignment
  void operator=(const pseudoECGElectro &);

public:
  //- Runtime type information
  TypeName("pseudoECGElectro");

  // Constructor
  pseudoECGElectro(const volScalarField &Vm, const dictionary &dict,
                   const volTensorField *conductivityPtr);

  // Destructor
  virtual ~pseudoECGElectro() = default;

  // Member Functions

  //- Compute Gima-Rudy pseudo-ECG for the current timestep
  virtual void compute();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
