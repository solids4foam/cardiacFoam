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
    SingleCellSolver

Description
    Single-cell electrophysiology solver based on OpenFOAM.

    This solver advances a single cardiac cell model in time using a
    run-time selectable ionic model. In contrast to electroActivationFoam,
    no spatial voltage PDE is solved; instead, the transmembrane voltage
    Vm is integrated directly as part of the ionic model ODE system.

    The solver is primarily intended for:
      - single-cell action potential simulations,
      - testing and validation of ionic models,
      - debugging and development of electrophysiology model components.

    External stimulation is prescribed via a stimulus protocol dictionary.
    Output is written as time series of states, algebraic variables and rates
    using ionicModelIO.

Author
    Simao Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    SingleCellSolver.C
\*---------------------------------------------------------------------------*/

#ifndef SingleCellSolver_H
#define SingleCellSolver_H

#include "electroModel.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"
#include "scalarField.H"
#include "overrideTypeName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electroModels
{

/*---------------------------------------------------------------------------*\
                        Class SingleCellSolver Declaration
\*---------------------------------------------------------------------------*/

class SingleCellSolver
:
    public electroModel
{
    //- Run-time selectable ionic model
    autoPtr<ionicModel> ionicModelPtr_;

    //- Optional ionic-model pre/post processor selected from the ionic model
    //  type
    autoPtr<electroVerificationModel> verificationModelPtr_;

    //- Processor-declared fields used for pre-processing data exchange
    wordList preProcessFieldNames_;
    PtrList<volScalarField> preProcessFields_;

    //- Processor-declared scratch fields used only for post-processing
    wordList postProcessFieldNames_;
    PtrList<volScalarField> postProcessFields_;

    //- Output file pointer
    autoPtr<OFstream> outputPtr_;

    //- Transmembrane voltage field
    volScalarField Vm_;

    //- Ionic model output fields used by lifecycle hooks
    PtrList<volScalarField> outFields_;

    //- Reused scratch field for the ionic current output channel
    scalarField dummyIonicCurrentField_;

    //- Disallow default bitwise copy construct
    SingleCellSolver(const SingleCellSolver&);

    //- Disallow default bitwise assignment
    void operator=(const SingleCellSolver&);


public:

    //- Runtime type information
    OverrideTypeName("singleCellSolver");

    // Constructors

        //- Construct from components
        SingleCellSolver
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    //- Destructor
    virtual ~SingleCellSolver() {}

    // Member Functions

        //- SingleCellSolver is ODE-only; deltaT is driven by controlDict.
        //  Override to prevent the base class from dereferencing the empty
        //  domainSystem_.
        virtual void setDeltaT(Time& runTime) override
        {
            physicsModel::setDeltaT(runTime);
        }

        //- Evolve the electro model
        virtual bool evolve() override;

        //- This function is called at the end of the simulation
        virtual void end() override;

        //- Return transmembrane voltage (needed by composite models)
        virtual const volScalarField& Vm() const override
        {
            return Vm_;
        }

        //- Return optional coupling signal provider (e.g. ionic model)
        virtual const ElectromechanicalSignalProvider* provider() const override
        {
            return ionicModelPtr_.operator->();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
