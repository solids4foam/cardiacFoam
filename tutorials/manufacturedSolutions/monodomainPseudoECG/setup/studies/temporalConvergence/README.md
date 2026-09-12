# temporalConvergence - monodomainPseudoECG

## Purpose

This study validates temporal convergence (timestep refinement) using the monodomainPseudoECG exact solution.

## Execution

From the repository root:

```bash
driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/temporalConvergence/sweep_temporal_convergence.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-temporal
```

The 12 cases are fixed-mesh timestep ladders: four levels at `N=640` for 1D
and 2D and four at `N=160` for 3D.  The pseudo-ECG verifier remains active so
field and functional temporal responses are archived together.  RKF45 uses
explicit baseline controls `absTol=1e-10` and `relTol=1e-8`; the companion
`odeToleranceControl` study repeats the finest level of each dimensional
ladder with tighter values.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
