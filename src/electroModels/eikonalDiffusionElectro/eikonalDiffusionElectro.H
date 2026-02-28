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
    eikonalDiffusionElectro

Description
    Solver for computing cardiac activation times by solving an anisotropic
    eikonal-type equation.

    The solver computes the activation time field psi(x), whose isocontours
    represent the propagation of the electrical depolarisation wavefront
    through cardiac tissue. The governing equation is a steady nonlinear
    eikonal–diffusion formulation derived as a simplification of the cardiac
    bidomain model, accounting for tissue anisotropy through a conductivity
    tensor.

    In its general form, the equation solved is:

        - div(M & grad(psi))
        + c0*sqrt(grad(psi) & M & grad(psi))
       == 1

    where:
        psi   is the activation time,
        c0    is a scaled conduction parameter (units T^{-1/2}), calibrated
              such that c0*sqrt(M) equals the prescribed physical fibre
              conduction velocity for a planar wavefront
        M     is the anisotropic conductivity (or metric) tensor constructed
              from myocardial fibre directions.

    The equation is solved iteratively using a finite-volume discretisation.
    Two solution strategies are provided: a simple Picard iteration in which
    the nonlinear eikonal metric term is treated explicitly, and an optional
    stabilised Picard (defect-correction) formulation in which a linearised
    advection operator, derived from the nonlinear eikonal term, is included
    implicitly to improve convergence. The selection between these approaches
    is controlled by the eikonalAdvectionDiffusionApproach switch; in both
    cases, the original steady eikonal–diffusion equation is recovered at
    convergence.

    Prescribed activation times may be imposed on selected regions or patches
    to represent electrical stimuli, while zero-gradient conditions are applied
    elsewhere.

    The resulting activation time field can be used directly for analysis of
    electrical wave propagation or as input to subsequent electrophysiology
    or electrocardiographic models.

    References
    Quarteroni, A., Lassila, T., Rossi, S., & Ruiz-Baier, R. (2017).
    Mathematical modelling of the human cardiovascular system: Data, numerical
    approximation, clinical applications. Springer, MS&A – Modeling, Simulation
    and Applications, Vol. 33.

Author
    Philip Cardiff, UCD.

SourceFiles
    eikonalDiffusionElectro.C

\*---------------------------------------------------------------------------*/

#ifndef eikonalDiffusionElectro_H
#define eikonalDiffusionElectro_H

#include "electroModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

/*---------------------------------------------------------------------------*\
                        Class eikonalDiffusionElectro Declaration
\*---------------------------------------------------------------------------*/

class eikonalDiffusionElectro
:
    public electroModel
{
    // Private data

        //- Activation time
        volScalarField psi_;

        // Gradient of activation time
        volVectorField gradPsi_;

        //- Indices of the cells to stimulate
        labelList stimulusCellIDs_;

        //- Surface-to-volume ratio
        const dimensionedScalar chi_;

        //- Membrane capacitance
        const dimensionedScalar Cm_;

        //- C-nductiivty tensor field
        volTensorField conductivity_;

        //- M tensor
        const volTensorField M_;

        //- M & gradPsi term
        volVectorField w_;

        //- Nonlinear term
        volScalarField G_;

        //- Reference "velocity"
        const dimensionedScalar c0_;

        //- Term related to the linearisation of c0*G
        volVectorField a_;

        //- Term related to the linearisation of c0*G
        volVectorField u_;

        //- phiU term
        surfaceScalarField phiU_;

        //- divPhiU term
        volScalarField divPhiU_;

        //- Solution algorithm flag
        //  If enabled, then an eikonal-advection-diffusion equation is solved,
        //  otherwise an eikonal-diffusion equation is solved. The eikonal-
        //  advection-diffusion approach should provide fasater convergence.
        const Switch eikonalAdvectionDiffusionApproach_;

    // Private Member Functions

        //- Initialise the conductivity
        tmp<volTensorField> initialiseConductivity() const;

        //- Disallow default bitwise copy construct
        eikonalDiffusionElectro(const eikonalDiffusionElectro&);

        //- Disallow default bitwise assignment
        void operator=(const eikonalDiffusionElectro&);


public:

    //- Runtime type information
    TypeName("eikonalDiffusionElectro");

    // Constructors

        //- Construct from components
        eikonalDiffusionElectro
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~eikonalDiffusionElectro()
        {}


    // Member Functions

        // Access

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime)
            {
                InfoInFunction
                    << "Setting deltaT and endTime to 1.0" << endl;
                runTime.setDeltaT(1.0);
                runTime.setEndTime(runTime.deltaT());
            }

            //- Evolve the electro model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
