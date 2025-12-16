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

Solver
    eikonalElectroActivationFoam

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

Authors
   Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // Create pimple loop object
    pimpleControl pimple(mesh);

    Info<< nl << "Solving eikonal-diffusion equation for psi" << nl << endl;

    // Set the time-step to 1 and increment the time
    runTime.setDeltaT(1.0);
    runTime++;

    while (pimple.loop())
    {
        // Update the gradient
        gradPsi = fvc::grad(psi);

        // Update w term
        w = M & gradPsi;

        // Update the nonlinear term
        G = sqrt((gradPsi & w) + smallG);

        // Solve eikonal equation
        if (eikonalAdvectionDiffusionApproach)
        {
            // Update a term
            a = w/G;

            // Update the u term
            u = c0*a;

            // Update the phiU term
            phiU = (fvc::interpolate(u) & mesh.Sf());

            // Construct the eikonal–diffusion equation (Eq. 6.14 in Quarteroni
            // et al.)
            //
            // The equation is written in a stabilised Picard (defect-
            // correction) form.
            // The implicit advection operator on the LHS is derived from a
            // linearisation of the nonlinear eikonal term c0*G, but is used
            // here purely to improve convergence. The corresponding explicit
            // advection term on the RHS ensures that, at convergence of the
            // outer iterations, the original eikonal–diffusion equation is
            // recovered exactly.
            //
            // The SuSp formulation is used to avoid diagonal weakening arising from
            // div(phiU), improving robustness of the linear solve.
            fvScalarMatrix psiEqn
            (
              - fvm::laplacian(M, psi)
              + fvm::div(phiU, psi) + fvm::SuSp(-fvc::div(phiU), psi)
             == one
              + fvc::div(phiU, psi) - fvc::div(phiU)*psi
              - c0*G
            );

            // Enforce psi to be zero for the stimulus cells
            psiEqn.setValues(stimulusCellIDs, 0.0);

            // Solve the system for psi
            psiEqn.solve("asymmetric_" + psi.name());
        }
        else
        {
            // Construct the eikonal-diffusion equation, where the nonlinear
            // term c0*G is treated entirely as an explicit deferred correction
            fvScalarMatrix psiEqn
            (
              - fvm::laplacian(M, psi)
              + c0*G
             == one
            );

            // Optional equation relaxation
            psiEqn.relax();

            // Enforce psi to be zero for the stimulus cells
            psiEqn.setValues(stimulusCellIDs, 0.0);

            // Solve the system for psi
            psiEqn.solve();
        }
    }

    Info<< nl;

    runTime.printExecutionTime(Info);

    Info<< "Writing " << psi.name() << " to " << runTime.value() << endl;

    psi.write();

    Info<< "End" << nl << endl;

    return 0;
}
