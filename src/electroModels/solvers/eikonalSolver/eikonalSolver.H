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
    EikonalSolver

Description
    Solver for computing cardiac activation times by solving an anisotropic
    eikonal-type equation.

    The solver computes the activation time field psi(x), whose isocontours
    represent the propagation of the electrical depolarisation wavefront
    through cardiac tissue. The governing equation is a steady nonlinear
    eikonal-diffusion formulation derived as a simplification of the cardiac
    bidomain model, accounting for tissue anisotropy through a conductivity
    tensor.

Author
    Philip Cardiff, UCD.

SourceFiles
    EikonalSolver.C
\*---------------------------------------------------------------------------*/

#ifndef EikonalSolver_H
#define EikonalSolver_H

#include "electroModel.H"
#include "overrideTypeName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electroModels
{

/*---------------------------------------------------------------------------*\
                        Class EikonalSolver Declaration
\*---------------------------------------------------------------------------*/

class EikonalSolver
:
    public electroModel
{
    //- Activation time
    volScalarField psi_;

    //- Transmembrane voltage (not solved here, but needed by the interface)
    volScalarField Vm_;

    //- Gradient of activation time
    volVectorField gradPsi_;

    //- Indices of the cells to stimulate
    labelList stimulusCellIDs_;

    //- Surface-to-volume ratio
    const dimensionedScalar chi_;

    //- Membrane capacitance
    const dimensionedScalar Cm_;

    //- Conductivity tensor field
    volTensorField conductivity_;

    //- M tensor
    const volTensorField M_;

    //- M & gradPsi term
    volVectorField w_;

    //- Nonlinear term
    volScalarField G_;

    //- Reference velocity
    const dimensionedScalar c0_;

    //- Terms related to the linearisation of c0*G
    volVectorField a_;
    volVectorField u_;

    //- phiU term
    surfaceScalarField phiU_;

    //- divPhiU term
    volScalarField divPhiU_;

    //- Solution algorithm flag
    const Switch eikonalAdvectionDiffusionApproach_;

    //- Initialise the conductivity
    tmp<volTensorField> initialiseConductivity() const;

    //- Disallow default bitwise copy construct
    EikonalSolver(const EikonalSolver&);

    //- Disallow default bitwise assignment
    void operator=(const EikonalSolver&);


public:

    //- Runtime type information
    OverrideTypeName("EikonalSolver");

    // Constructors

        //- Construct from components
        EikonalSolver
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    //- Destructor
    virtual ~EikonalSolver() {}

    // Member Functions

        //- Update the size of the time-step
        virtual void setDeltaT(Time& runTime) override
        {
            InfoInFunction << "Setting deltaT and endTime to 1.0" << endl;
            runTime.setDeltaT(1.0);
            runTime.setEndTime(runTime.deltaT());
        }

        //- Evolve the electro model
        virtual bool evolve() override;

        //- Return transmembrane voltage (needed by composite models)
        virtual const volScalarField& Vm() const override
        {
            return Vm_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
