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
    monoDomainElectro

Description
    Generic model-agnostic solver for cardiac electrophysiology based on the
    monodomain reaction-diffusion equation.

Author
    Simao Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    monoDomainElectro.C
\*---------------------------------------------------------------------------*/

#ifndef monoDomainElectro_H
#define monoDomainElectro_H

#include "electroModel.H"
#include "ionicModel.H"
#include "ecgModel.H"
#include "purkinjeModel.H"
#include "electroModelsPrePostProcessor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electroModels
{

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


/*---------------------------------------------------------------------------*\
                        Class monoDomainElectro Declaration
\*---------------------------------------------------------------------------*/

class monoDomainElectro
:
    public electroModel
{
    //- Transmembrane voltage
    volScalarField Vm_;

    //- Ionic current
    volScalarField Iion_;

    //- External stimulus current
    volScalarField externalStimulusCurrent_;

    //- List of external stimuli
    List<stimulus> externalStimulus_;

    //- Cardiac properties dictionary
    const dictionary& cardiacProperties_;

    //- Surface-to-volume ratio
    const dimensionedScalar chi_;

    //- Membrane capacitance
    const dimensionedScalar Cm_;

    //- Conductivity tensor field
    volTensorField conductivity_;

    //- Activation time field
    volScalarField activationTime_;

    //- Per-cell flag indicating if the activation time is still unset
    boolList calculateActivationTime_;

    //- Ionic model output fields
    PtrList<volScalarField> outFields_;

    //- Run-time selectable ionic model
    autoPtr<ionicModel> ionicModelPtr_;

    //- Optional ionic-model pre/post processor
    autoPtr<electroModelsPrePostProcessor> prePostProcessor_;

    //- Optional ECG model selected from the nested ECG configuration
    autoPtr<ecgModel> ecgModelPtr_;

    //- Processor-declared fields used for pre-processing data exchange
    wordList preProcessFieldNames_;
    PtrList<volScalarField> preProcessFields_;

<<<<<<< HEAD
    //- Flag for tracking if the time-step size has been set
    bool setDeltaT_;
=======
  //- Optional 1D Purkinje Network model
  autoPtr<purkinjeModel> purkinjeModelPtr_;

  //- Processor-declared fields used for pre-processing data exchange
  wordList preProcessFieldNames_;
  PtrList<volScalarField> preProcessFields_;
>>>>>>> eb1d0e9 (feat: wire purkinjeModelPtr_ into monoDomainElectro)

    //- Frequency for writing time-step information to standard output
    label infoFrequency_;

    //- Initialise the conductivity
    tmp<volTensorField> initialiseConductivity() const;

    //- Update the external stimulus current field
    void updateExternalStimulusCurrent
    (
        volScalarField& externalStimulusCurrent,
        const List<stimulus>& externalStimulus,
        const scalar t0
    ) const;

    //- Update the activation time
    void updateActivationTime
    (
        volScalarField& activationTime,
        boolList& calculateActivationTime,
        const volScalarField& Vm
    ) const;

    //- Evolve the electro model using an explicit algorithm
    virtual bool evolveExplicit();

    //- Evolve the electro model using an implicit algorithm
    virtual bool evolveImplicit();

    //- Disallow default bitwise copy construct
    monoDomainElectro(const monoDomainElectro&);

    //- Disallow default bitwise assignment
    void operator=(const monoDomainElectro&);


protected:

    //- Protected constructor for subclasses
    monoDomainElectro
    (
        const word& type,
        Time& runTime,
        const word& region = dynamicFvMesh::defaultRegion
    );

    //- Return const access to the conductivity tensor
    const volTensorField& conductivity() const
    {
        return conductivity_;
    }

    //- Re-read Vm from the current time directory
    void readVm();


public:

    //- Initialise a conductivity tensor field
    static tmp<volTensorField> initialiseConductivityTensor
    (
        const word& fieldName,
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Runtime type information
    TypeName("monoDomainElectro");

    // Constructors

        //- Construct from components
        monoDomainElectro
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    //- Destructor
    virtual ~monoDomainElectro() {}

    // Member Functions

        //- Update the size of the time-step
        virtual void setDeltaT(Time& runTime);

        //- Evolve the electro model
        virtual bool evolve();

        //- Read electroProperties dictionary
        virtual bool read();

        //- Return transmembrane voltage (needed by composite models)
        virtual const volScalarField& Vm() const
        {
            return Vm_;
        }

        //- Return optional coupling signal provider (e.g. ionic model)
        virtual const CouplingSignalProvider* provider() const
        {
            return ionicModelPtr_.operator->();
        }

        //- Finalise electro and optional ECG reporting
        virtual void end();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
