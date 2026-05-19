# tutorials architecture

This folder contains reference and regression cases for `cardiacFoam`.

## Current tutorial cases

- `singleCellprotocols/singleCell` : single-point ODE workflow (`singleCellSolver`)
- `singleCellprotocols/ionicHeterogeneityProbe` : meshless Bueno-Orovio
  transmural heterogeneity probe with 2D/3D plotting
- `singleCellprotocols/restitutionCurves_s1s2Protocol` : S1-S2 pacing sweeps
  (`singleCellSolver`)
- `NiedererEtAl2011/NiedererEtAl2011verification` : slab verification workflow
  (`myocardiumSolver monodomainSolver`)
- `NiedererEtAl2011/monodomainPurkinjeNiedererEtAl2011` : Niederer slab with a
  small 1D Purkinje network
- `NiedererEtAl2011/electroMechanicalNiedererEtAl2011` : electromechanical
  Niederer slab using `electroMechanicalModel` with monodomain electrophysiology
- `manufacturedSolutions/monodomainPseudoECG` : spatial manufactured-solution
  verification with pseudo-ECG (`monodomainFDAManufactured`)
- `manufacturedSolutions/bidomain` : spatial manufactured-solution verification
  (`bidomainFDAManufactured`)
- `manufacturedSolutions/bathBidomain` : bidomain-with-bath manufactured
  verification (`bathBidomainFDAManufactured`, `torsoECG`)

## Common script pattern

Most runnable cases provide:

- `Allrun` : run simulation and sometimes post-process
- `Allclean` : remove generated output
- optional `regressionTest.sh` or `runRegressionTest.sh` : case-local numerical
  checks

## Cross-case regression entrypoint

From `tutorials/`:

```bash
./Alltest-regression
```

Each covered tutorial owns its regression script and reference file locally;
there is no separate shared regression tree.

## Python automation integration

The shared driver (`foamctl` / `openfoam_driver`) maps these case folders to
tutorial specs and can run parameter sweeps with reproducible manifests and
post-processing.
