# coupledConvergence - monodomain1D3D

## Purpose

This study validates spatial and/or temporal convergence for the coupled 1D-3D monodomain1D3D solver.

## Execution

From the repository root, use the current wrapper and one of the study specs:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/monodomain1D3D/setup/studies/coupledConvergence/sweep_active.json \
    --output-dir .tmp/driverfoam/monodomain1D3D-active
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/monodomain1D3D/setup/studies/coupledConvergence/sweep_active.json \
    --output-dir .tmp/driverfoam/monodomain1D3D-active
```

The `bidirectional` and `decoupled` specs currently fail before OpenFOAM because
their overrides target `couplingMode` and `rPvj` without the current override
scope. Those keys live under
`monodomainSolverCoeffs.domainCouplings.couplingA`; fix that
spec/dictionary contract before running those variants.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
