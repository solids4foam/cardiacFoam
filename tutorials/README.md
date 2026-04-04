# tutorials architecture

This folder contains reference and regression cases for `cardiacFoam`.

## Current tutorial cases

- `singleCell` : single-point ODE workflow (`singleCell`)
- `NiedererEtAl2012` : slab verification workflow (`monodomain`)
- `ECG` : monodomain + ECG output workflow (`electroModel` with nested `ECG`)
- `manufacturedFDA` : manufactured-solution verification (`tmanufacturedFDA`)
- `restitutionCurves_s1s2Protocol` : S1-S2 pacing sweeps (single-cell)
- `vortexDynamics` : 2D wave dynamics (monodomain)

## Common script pattern

Most cases provide:

- `Allrun` : run simulation (and sometimes post-process)
- `Allclean` : remove generated output
- optional `runRegressionTest.sh` : case-specific numerical checks

## Cross-case regression entrypoint

From `tutorials/`:

```bash
./Alltest-regression
```

Defaults to a smoke set and supports:

- `smoke`
- `all`
- `<case>`
- `<case>_parallel`

## Python automation integration

The shared driver (`foamctl` / `openfoam_driver`) maps these case folders to tutorial
specs and can run parameter sweeps with reproducible manifests and post-processing.
