# temporalConvergence - monodomainPseudoECG

## Purpose

This study validates temporal convergence (timestep refinement) using the monodomainPseudoECG exact solution.

## Execution

From the repository root:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/temporalConvergence/sweep_temporal_convergence.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-temporal
```

This spec currently fails during materialisation for the same temporal
factory-cardinality issue as the bidomain temporal study. Repair the JSON or
factory contract before using `sweep-run`.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
