# tutorials architecture

This folder contains reference and regression cases for `cardiacFoam`.

## Current tutorial cases

- `singleCell` : single-point ODE workflow (`SingleCellSolver`)
- `Niederer/` : grouped slab verification and Purkinje coupling cases
  - `NiedererEtAl2012verification` : slab verification workflow (`MonoDomainSolver`)
  - `NiedererEtAl2012MonodomainPurkinje` : Niederer slab extended with a small 1D Purkinje network
  - `NiedererEtAl2012EikonalPurkinje` : Niederer slab configured for 3D eikonal + 1D Purkinje activation tests
- `ECG` : monodomain + ECG output workflow (`electroModel` with nested `ECG`)
- `manufacturedSolutions/` : grouped manufactured-solution verification cases
  - `monodomainPseudoECG` : spatial manufactured-solution verification with pseudo-ECG (`monodomainFDAManufactured`)
  - `bidomain` : spatial manufactured-solution verification (`bidomainFDAManufactured`)
  - `singleCellMonodomain` : single-cell manufactured-solution verification
  - `singleCellBidomain` : single-cell manufactured-solution verification
- `regressionTests/` : shared regression assets separated from runnable tutorial folders
  - `NiedererEtAl2012` : Niederer slab regression overrides and checkpoints
  - `singleCell` : single-cell regression overrides and checkpoints
- `restitutionCurves_s1s2Protocol` : S1-S2 pacing sweeps (single-cell)
- `vortexDynamics` : 2D wave dynamics (monodomain)

## Common script pattern

Most runnable cases provide:

- `Allrun` : run simulation (and sometimes post-process)
- `Allclean` : remove generated output
- optional `runRegressionTest.sh` : thin wrapper to case-specific numerical checks

## Cross-case regression entrypoint

From `tutorials/`:

```bash
./regressionTests/Alltest-regression
```

Defaults to a smoke set and supports:

- `smoke`
- `all`
- `<case>`
- `<case>_parallel`

The regression inputs themselves live under `tutorials/regressionTests/`, while
the runnable tutorials keep local wrappers for convenience.

## Python automation integration

The shared driver (`foamctl` / `openfoam_driver`) maps these case folders to tutorial
specs and can run parameter sweeps with reproducible manifests and post-processing.
