# tetConvergence - bidomain

## Purpose

This study validates spatial convergence on unstructured tetrahedral meshes using the bidomain exact solution.

## Execution

From the repository root, use the current wrapper:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/tetConvergence/sweep_tet_generic.json \
    --output-dir .tmp/driverfoam/bidomain-tet
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/tetConvergence/sweep_tet_generic.json \
    --output-dir .tmp/driverfoam/bidomain-tet
```

Resolve the runtime preflight before expecting OpenFOAM execution.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
