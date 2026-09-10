# nonlinearControl - eikonalECG

## Purpose

This twelve-case control reruns the N = 10, 20, 40, and 80 least-squares,
generic-Delaunay tetrahedral cases for all three headline eikonal regimes:
axis-aligned conductivity, rotated conductivity without the advection-diffusion
approach, and rotated conductivity with it. The eikonal outer fixed-point
stopping tolerance is tightened from 1e-8 to 1e-10. It isolates the nonlinear
stopping criterion while retaining the same manufactured solution, conductivity,
mesh family, gradient scheme, and ECG settings as the matching rows of
../tetConvergence/sweep_tet_generic.json.

The Paper I release matrix deliberately uses one unstructured mesh-generation
family: the generic Gmsh Delaunay mesh defined in
../tetConvergence/box.geo.template. Frontal/Netgen is not part of this control
or the planned release reruns.

## Execution

First materialize and inspect the twelve driverFOAM cases:

    driverFoam sweep-plan \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/nonlinearControl/sweep_tet_outer_tolerance.json

Then execute the same manifest with driverFOAM:

    driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/nonlinearControl/sweep_tet_outer_tolerance.json

Compare the archived activation-time and ECG metrics against the matching rows
of the generic tetrahedral baseline, and retain the complete four-level
sequence to test whether the conclusion holds over the fitted convergence
window. Compute each pairwise order using the measured cell-count scale
\(h=n_{\mathrm{cells}}^{-1/3}\), rather than assuming that the nominal Gmsh size
parameter halves exactly. The control supports the claim that the reported
spatial error is not limited by the eikonal outer stopping criterion.

## Outputs

Generated case outputs and archives are written under this study's local
results/ directory and must not be committed.
