# cartesianConvergence - eikonalECG

## Purpose

This study validates spatial convergence on structured hexahedral (cartesian) meshes using the eikonalECG exact solution.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json
python3 applications/scripts/paperI_results/aggregate.py eikonal_cartesian
```

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
