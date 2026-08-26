# temporalConvergence - bidomain

## Purpose

This study validates temporal convergence (timestep refinement) using the bidomain exact solution.

## Execution

From the repository root:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/temporalConvergence/sweep_temporal_convergence.json \
    --output-dir .tmp/driverfoam/bidomain-temporal
```

This spec currently fails during materialisation because its temporal `dt`
ladder is interpreted as four factory cases for each dimension, while the
entry-based sweep contract requires one case per expanded axis combination.
Repair that JSON/factory contract before using `sweep-run`.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
