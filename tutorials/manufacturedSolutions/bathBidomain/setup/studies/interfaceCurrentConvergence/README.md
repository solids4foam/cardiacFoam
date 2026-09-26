# interfaceCurrentConvergence - bathBidomain

## Purpose

Sweeps assembled interface-current rows for `@tbl-bath-bidomain-tet`: the conformal tetrahedral bath-bidomain system at `N=10,20,40,80`, `matchedSubmesh` flux assembly, comparing the `unweightedHarmonic` and `distanceWeightedHarmonic` interface-conductivity interpolation methods. Both write the assembled-current metrics computed by the `bathBidomainInterfaceMetrics` function object (invoked `-latestTime` by driverFOAM's own tet workflow DAG) with `snGrad corrected` (the case's committed default).

Replaces the former `setup/studies/tetConvergence/run_parallel_interface_sweep.sh` bash script — mesh generation, `checkMesh`, and the tet electroProperties/fvSchemes overlay swap are now handled directly by driverFOAM's own `manufacturedBathBidomain` tet workflow DAG, not by hand-rolled bash.

## Execution

```bash
driverFoam sweep-run --spec setup/studies/interfaceCurrentConvergence/sweep_tet_unweightedHarmonic.json
driverFoam sweep-run --spec setup/studies/interfaceCurrentConvergence/sweep_tet_distanceWeightedHarmonic.json
```

Each method is its own spec (rather than a single sweep with a `method` axis) because the sweep engine's `caseId`/`output_dir_name` templating can only reference scalar/list independent axes; `electro_property_overrides` (the mechanism that sets `interfaceConductivityInterpolation`) is a per-case dict and isn't safe to reference there, so it's fixed in each spec's `base` instead.

## Status

Verified with a real `driverFoam sweep-run` at `N=10` (`unweightedHarmonic`, 2026-08-19): the case meshes, solves, and writes `bathBidomainInterfaceMetrics.csv` under the sweep's own archive layout (`<case_root>/<caseId>/setup/studies/interfaceCurrentConvergence/results/sweepCases...`). `N=20,40,80` and the `distanceWeightedHarmonic` spec use the same mechanism and are not yet run. `bath_predictor_corrector: true` in `base` enables predictor-corrector coupling.

## Tracking & Outputs

All generated outputs land under the local `results/` folder, gitignored. Do not commit generated OpenFOAM data.
