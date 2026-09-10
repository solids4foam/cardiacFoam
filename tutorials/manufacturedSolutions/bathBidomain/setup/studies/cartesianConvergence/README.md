# cartesianConvergence - bathBidomain

## Purpose

This study validates spatial convergence on structured hexahedral (cartesian) meshes using the bathBidomain exact solution.

## Execution

From the repository root, use the current wrapper:

```bash
driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/bathBidomain/setup/studies/cartesianConvergence/sweep_hex_convergence.json \
    --output-dir .tmp/driverfoam/bathBidomain-cartesian
driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bathBidomain/setup/studies/cartesianConvergence/sweep_hex_convergence.json \
    --output-dir .tmp/driverfoam/bathBidomain-cartesian
```

Resolve the runtime preflight before expecting OpenFOAM execution.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
