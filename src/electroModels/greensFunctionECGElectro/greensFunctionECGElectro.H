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

Class
    greensFunctionECGElectro

Description
    Extends monoDomainElectro to compute body-surface ECG potentials via
    the Green's function (volume integral) approach.

    After each monodomain time step the current source density is computed:

        Is(x,t) = -div( Gi . grad(Vm) )

    and the potential at each electrode position P is evaluated as:

        phi(P,t) = 1/(4*pi*sigmaT) * sum_c [ Is_c * V_c / |C_c - P| ]

    where the sum runs over all heart mesh cells and is made parallel-safe
    via Foam::reduce.

    Gi is initialised by first attempting to read a volTensorField named
    "Gi" from disk (READ_IF_PRESENT). If not found, a uniform tensor is
    built from the ECG.Gi symmTensor entry in the dict, following the same
    pattern as monoDomainElectro::initialiseConductivity().

    All ECG parameters are read from an ECG subdict inside the model's
    <type>Coeffs dictionary in electroProperties:

        greensFunctionECGElectroCoeffs
        {
            // --- standard monodomain parameters ---
            chi   ...;
            Cm    ...;
            conductivity ...;
            ionicModel ...;
            // ...

            ECG
            {
                // symmTensor (6 components: xx xy xz yy yz zz)
                // OR a volTensorField named "Gi" present in the time directory
                Gi      [-1 -3 3 0 0 2 0] (0.17 0 0  0.019 0  0.019);

                sigmaT  [-1 -3 3 0 0 2 0] 0.2;

                electrodes
                {
                    LA  ( 0.050  0.120  0.060);
                    RA  (-0.050  0.120  0.060);
                    LL  ( 0.000 -0.150  0.000);
                }
            }
        }

    Output: postProcessing/ECG.dat — one row per time step:
        # time  LA  RA  LL ...

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    greensFunctionECGElectro.C

\*---------------------------------------------------------------------------*/

#ifndef greensFunctionECGElectro_H
#define greensFunctionECGElectro_H

#include "monoDomainElectro.H"
#include "OFstream.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

/*---------------------------------------------------------------------------*\
                   Class greensFunctionECGElectro Declaration
\*---------------------------------------------------------------------------*/

class greensFunctionECGElectro
:
    public monoDomainElectro
{
    // Private data

        //- Current source density  Is = -div(Gi . grad(Vm))  [A/m^3]
        volScalarField Is_;

        //- Intracellular conductivity tensor  [S/m]
        volTensorField Gi_;

        //- Torso conductivity (uniform, isotropic)  [S/m]
        const dimensionedScalar sigmaT_;

        //- Electrode names (empty if not configured)
        wordList electrodeNames_;

        //- Electrode positions (empty if not configured)
        List<vector> electrodePositions_;

        //- Output stream for ECG.dat (null if electrodes not configured)
        autoPtr<OFstream> outputPtr_;

        //- Torso surface for body-surface potential map (null if not configured)
        autoPtr<triSurface> torsoSurfacePtr_;

        //- Face centres of the torso surface
        pointField torsoFaceCentres_;

        //- Output directory for VTK surface files
        fileName torsoVtkDir_;


    // Private Member Functions

        //- Initialise Gi: read from disk if present, else use dict constant
        tmp<volTensorField> initialiseGi() const;

        //- Read electrode positions from the ECG subdict
        void readElectrodes(const dictionary& ecgDict);

        //- Write column header to ECG.dat
        void writeHeader();

        //- Load torso surface STL if "torsoSurface" key is present in ecgDict
        void loadTorsoSurface(const dictionary& ecgDict);

        //- Write torso surface potential as VTK ASCII POLYDATA
        void writeTorsoVtk(const scalarList& phiTorso) const;

        //- Disallow default bitwise copy construct
        greensFunctionECGElectro(const greensFunctionECGElectro&);

        //- Disallow default bitwise assignment
        void operator=(const greensFunctionECGElectro&);


public:

    //- Runtime type information
    TypeName("greensFunctionECGElectro");


    // Constructors

        //- Construct from components
        greensFunctionECGElectro
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~greensFunctionECGElectro()
        {}


    // Member Functions

        //- Evolve: monodomain step followed by ECG integral evaluation
        virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
