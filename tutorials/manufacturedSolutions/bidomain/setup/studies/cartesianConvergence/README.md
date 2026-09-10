# cartesianConvergence - bidomain

## Purpose

This study validates spatial convergence on structured hexahedral (cartesian) meshes using the bidomain exact solution.

## Execution

From the repository root, plan and run the sweep with:

```bash
driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/cartesianConvergence/sweep_hex_convergence.json \
    --output-dir .tmp/driverfoam/bidomain-cartesian
driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/cartesianConvergence/sweep_hex_convergence.json \
    --output-dir .tmp/driverfoam/bidomain-cartesian
```

Resolve the runtime preflight before expecting OpenFOAM execution; see
`tutorials/manufacturedSolutions/README.md`.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
