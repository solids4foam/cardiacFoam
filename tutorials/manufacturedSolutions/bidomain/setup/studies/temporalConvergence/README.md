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
both 1D and 2D. The configured `sbdf2` coupling is archived for every case.

Corrected 2026-09-26: this study's spec used to also carry
`ode_abs_tolerance`/`ode_rel_tolerance` (`1e-10`/`1e-8`), presented above as
applied RKF45 baseline controls. Neither key was ever read by `make_spec`,
so they had no effect on any run; they are removed from the spec (owner,
plan §5g Q9). The case's own `constant/electroProperties` sets `solver
RKF45` with no tolerance overrides.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
