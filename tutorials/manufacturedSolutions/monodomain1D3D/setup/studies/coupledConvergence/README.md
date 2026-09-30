# coupledConvergence - monodomain1D3D

## Purpose

This study validates spatial and/or temporal convergence for the coupled 1D-3D monodomain1D3D solver.

## Execution

From the repository root, use the current wrapper and one of the study specs. Spec: `tutorials/manufacturedSolutions/monodomain1D3D/setup/studies/coupledConvergence/sweep_active.json`.

```bash
[omnidriver command to run]
```

The `bidirectional` and `decoupled` specs use fully scoped
`$ELECTRO_MODEL_COEFFS.domainCouplings.couplingA` overrides so they remain
valid if the active myocardium-solver coefficient dictionary is renamed.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
