# tetConvergence - bidomain

## Purpose
This study validates spatial convergence on unstructured tetrahedral meshes using the bidomain exact solution.

## Execution
Execute the sweep using `foamctl run --sweep sweep_tet_generic.json` or by invoking the local runner script.

## Tracking & Outputs
All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
