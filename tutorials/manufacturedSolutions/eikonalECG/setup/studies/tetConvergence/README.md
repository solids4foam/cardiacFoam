# tetConvergence - eikonalECG

## Purpose

This study validates spatial convergence on one reproducible unstructured
tetrahedral mesh family: the generic Gmsh Delaunay mesh defined in
box.geo.template. It contains the full N = 10, 20, 40, 80 baseline matrix for
the axis and rotated conductivities, the configured advection-diffusion
variants, and the Gauss-linear/least-squares comparisons used to interpret the
eikonalECG exact solution.

## Execution

First materialize the 16 baseline cases and inspect the generated plans:

    applications/scripts/driverFoam/bin/driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/tetConvergence/sweep_tet_generic.json

Then run the same manifest with driverFOAM:

    applications/scripts/driverFoam/bin/driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/tetConvergence/sweep_tet_generic.json

The previous reference to applications/scripts/paperI_results/aggregate.py was
removed because that script is not present in this repository. Use the
driverFOAM sweep manifest and its archived verification outputs as the canonical
run record; add an analysis command only alongside a versioned in-repository
postprocessor.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
