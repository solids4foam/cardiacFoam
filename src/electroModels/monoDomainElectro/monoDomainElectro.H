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
    monoDomainElectro

Description
    Generic model-agnostic solver for cardiac electrophysiology
    based on the monodomain reaction–diffusion equation. The solver:

    - Delegates all ionic-model state indexing and ODE evaluation to the
      run-time selectable ionicModel.

    - Supports multiple time-integration strategies (explicit, implicit).

    - Stores ionic state vectors externally (one N-state vector per cell),
      ensuring complete separation between solver logic and model detail.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    monoDomainElectro.C

\*---------------------------------------------------------------------------*/

#ifndef monoDomainElectro_H
#define monoDomainElectro_H

#include "electroModel.H"
#include "ionicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

// Stimulus helper struct for storing info related to external stimulii
struct stimulus
{
    //- Cell indices to be stimulated
    labelList cellIDs_;

    //- Start time for stimulation
    scalar startTime_;

    //- Duration of stimulation
    scalar duration_;

    //- Intensity
    scalar intensity_;

    stimulus()
    :
        cellIDs_(0),
        startTime_(0),
        duration_(0),
        intensity_(0)
    {}
};


/*---------------------------------------------------------------------------* \
                        Class monoDomainElectro Declaration
\*---------------------------------------------------------------------------*/

class monoDomainElectro
:
    public electroModel
{
    // Private data

        //- Transmembrane voltage
        volScalarField Vm_;

        //- Ionic currnet
        volScalarField Iion_;

        //- External stimulus current
        volScalarField externalStimulusCurrent_;

        //- List of external stimulii
        List<stimulus> externalStimulus_;

        //- Cardiac properties dictionary
        const dictionary& cardiacProperties_;

        //- Surface-to-volume ratio
        const dimensionedScalar chi_;

        //- Membrane capacitance
        const dimensionedScalar Cm_;

        //- C-nductiivty tensor field
        volTensorField conductivity_;

        //- Field to record the activation time, which is the time at which
        //  each cell reaches or exceeds a transmembrane voltage of 0.0
        volScalarField activationTime_;

        //- Field for tracking if the activation time has been calculated yet in
        //  a cell
        boolList calculateActivationTime_;

        //- Activation velocity
        volVectorField activationVelocity_;

        //- Ionic model output fields
        PtrList<volScalarField> outFields_;

        //- Run-time selectable ionic model
        autoPtr<ionicModel> ionicModelPtr_;


    // Private Member Functions

        //- Initialise the conductivity
        tmp<volTensorField> initialiseConductivity() const;

        //- Disallow default bitwise copy construct
        monoDomainElectro(const monoDomainElectro&);

        //- Disallow default bitwise assignment
        void operator=(const monoDomainElectro&);


public:

    //- Runtime type information
    TypeName("monoDomainElectro");

    // Constructors

        //- Construct from components
        monoDomainElectro
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~monoDomainElectro()
        {}


    // Member Functions

        // Edit

            //- Evolve the electro model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
