# temporalConvergence - bidomain

## Purpose

This study validates temporal convergence (timestep refinement) using the bidomain exact solution.

## Execution

From the repository root. Spec: `tutorials/manufacturedSolutions/bidomain/setup/studies/temporalConvergence/sweep_temporal_convergence.json`.

```bash
[omnidriver command to run]
```

The eight cases are fixed-mesh timestep ladders: four levels at `N=640` in
both 1D and 2D. The configured `sbdf2` coupling is archived for every case.

The case's own `constant/electroProperties` sets `solver RKF45` with no
tolerance overrides.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
