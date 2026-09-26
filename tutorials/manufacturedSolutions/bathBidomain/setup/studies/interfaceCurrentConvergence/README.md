# interfaceCurrentConvergence - bathBidomain

## Purpose

Sweeps assembled interface-current rows for `@tbl-bath-bidomain-tet`: the conformal tetrahedral bath-bidomain system at `N=10,20,40,80`, `matchedSubmesh` flux assembly, comparing the `unweightedHarmonic` and `distanceWeightedHarmonic` interface-conductivity interpolation methods. Both write the assembled-current metrics computed by the `bathBidomainInterfaceMetrics` function object (invoked `-latestTime` by the omniD `manufacturedBathBidomain` record's tet route) with `snGrad corrected` (the case's committed default).

Replaces the former `setup/studies/tetConvergence/run_parallel_interface_sweep.sh` bash script — mesh generation and `checkMesh` are the record's tet route, not hand-rolled bash.

## Execution

```bash
driverFoam sweep-run --spec setup/studies/interfaceCurrentConvergence/sweep_tet_unweightedHarmonic.json
driverFoam sweep-run --spec setup/studies/interfaceCurrentConvergence/sweep_tet_distanceWeightedHarmonic.json
```

Each method is its own spec. `sweep_tet_unweightedHarmonic.json` sets `bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation` in its `base`; `distanceWeightedHarmonic` is the case's own value, so its spec states nothing.

## Status

Verified with a real `driverFoam sweep-run` at `N=10` (`unweightedHarmonic`, 2026-08-19): the case meshes, solves, and writes `bathBidomainInterfaceMetrics.csv` under the sweep's own archive layout (`<case_root>/<caseId>/setup/studies/interfaceCurrentConvergence/results/sweepCases...`). `N=20,40,80` and the `distanceWeightedHarmonic` spec use the same mechanism and are not yet run. Predictor-corrector coupling is the case's own `bathPredictorCorrector yes`. (Corrected 2026-09-26: this said `bath_predictor_corrector: true` in `base`, a key the record study no longer carries, since it restated that value.)

## Tracking & Outputs

All generated outputs land under the local `results/` folder, gitignored. Do not commit generated OpenFOAM data.
