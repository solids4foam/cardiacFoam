# temporalConvergence - bidomain

## Purpose

This study validates temporal convergence (timestep refinement) using the bidomain exact solution.

## Execution

From the repository root:

```bash
driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/temporalConvergence/sweep_temporal_convergence.json \
    --output-dir .tmp/driverfoam/bidomain-temporal
```

The eight cases are fixed-mesh timestep ladders: four levels at `N=640` in
both 1D and 2D. The configured `sbdf2` coupling and adaptive RKF45 baseline
controls (`absTol=1e-10`, `relTol=1e-8`) are archived for every case. The
companion `odeToleranceControl` study repeats the finest level of each ladder
with tighter ODE tolerances.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
