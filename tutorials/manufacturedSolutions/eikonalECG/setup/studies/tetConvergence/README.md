# tetConvergence - eikonalECG

## Purpose

This study validates spatial convergence on unstructured tetrahedral meshes using the eikonalECG exact solution.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/tetConvergence/sweep_tet_generic.json
python3 applications/scripts/paperI_results/aggregate.py eikonal_tet_generic
```

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
