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

        phi_pseudo(P,t) = -sum_c [ (conductivity . grad(Vm))_c . r_vec * V_c / |r|^3 ]

    where r_vec = C_c - P and the sum is made parallel-safe via Foam::reduce.

    Dictionary (inside pseudoECGElectroCoeffs):

        pseudoECGElectroCoeffs
        {
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
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                     Class pseudoECGElectro Declaration
\*---------------------------------------------------------------------------*/

class pseudoECGElectro : public ecgModel {
  // Private data

  //- Output stream for pseudoECG.dat
  autoPtr<OFstream> outputPtr_;

  //- Optional manufactured-solution diagnostics stream
  autoPtr<OFstream> manufacturedOutputPtr_;

  //- Enable manufactured-reference comparisons after each ECG step
  Switch manufacturedEnabled_;

  //- Active manufactured dimension used by the meshless reference integral
  label manufacturedDimension_;

  //- High-order reference quadrature order
  label manufacturedReferenceQuadratureOrder_;

  //- Lower-order quadrature orders used for consistency checks
  List<label> manufacturedCheckQuadratureOrders_;

  //- Cached spatial manufactured reference for the higher-order quadrature
  List<scalar> manufacturedReferenceSpatialValues_;

  //- Cached spatial manufactured references for the lower-order quadratures
  List<List<scalar>> manufacturedCheckSpatialValues_;

  //- Conductivity tensor used to build the current manufactured cache
  tensor manufacturedCachedConductivity_;

  //- True when the cached manufactured spatial references are up to date
  bool manufacturedSpatialCacheValid_;

  //- Number of manufactured samples processed so far
  label manufacturedSampleCount_;

  //- Running norms for |numeric - reference quadrature|
  List<scalar> referenceErrorL1Sum_;
  List<scalar> referenceErrorL2Sum_;
  List<scalar> referenceErrorLinf_;

  //- Running norms for |numeric - check quadrature| for each qCheck
  List<List<scalar>> checkErrorL1Sum_;
  List<List<scalar>> checkErrorL2Sum_;
  List<List<scalar>> checkErrorLinf_;

  //- Running norms for |reference quadrature - check quadrature| for each qCheck
  List<List<scalar>> referenceDeltaL1Sum_;
  List<List<scalar>> referenceDeltaL2Sum_;
  List<List<scalar>> referenceDeltaLinf_;

  //- Avoid rewriting the manufactured summary file
  bool manufacturedSummaryWritten_;

  // Private Member Functions

  void readManufacturedConfig(const dictionary &dict);

  void initialiseManufacturedOutput();

  void resizeManufacturedCheckStorage();

  void invalidateManufacturedReferenceCache();

  void rebuildManufacturedReferenceCache(const tensor &referenceConductivity);

  void computePseudoECGValues(List<scalar> &values) const;

  void updateManufacturedStatistics(const List<scalar> &numericValues,
                                    const List<List<scalar>> &checkReferenceValues,
                                    const List<scalar> &referenceValues);

  void writeManufacturedSummary();

  //- Disallow copy construct
  pseudoECGElectro(const pseudoECGElectro &);

  //- Disallow assignment
  void operator=(const pseudoECGElectro &);

public:
  //- Runtime type information
  TypeName("pseudoECGElectro");

  // Constructor
  pseudoECGElectro(const volScalarField &Vm, const dictionary &dict,
                   const volTensorField &conductivity);

  // Destructor
  virtual ~pseudoECGElectro() = default;

  // Member Functions

  //- Compute Gima-Rudy pseudo-ECG for the current timestep
  virtual void compute();

  //- Refresh coefficients while preserving the existing output layout
  virtual bool read(const dictionary &dict);

  //- Flush manufactured diagnostics at the end of the run
  virtual void end();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
