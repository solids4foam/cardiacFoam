# temporalConvergence - bidomain

## Purpose

This study validates temporal convergence (timestep refinement) using the bidomain exact solution.

## Execution

Execute the sweep using `foamctl run --sweep sweep_temporal_convergence.json` or by invoking the local runner script.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
