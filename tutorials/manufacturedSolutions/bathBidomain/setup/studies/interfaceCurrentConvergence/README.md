# interfaceCurrentConvergence - bathBidomain

## Purpose

Sweeps assembled interface-current rows for `@tbl-bath-bidomain-tet`: the conformal tetrahedral bath-bidomain system at `N=10,20,40,80`, `matchedSubmesh` flux assembly, comparing the `unweightedHarmonic` and `distanceWeightedHarmonic` interface-conductivity interpolation methods. Both write the assembled-current metrics computed by the `bathBidomainInterfaceMetrics` function object (invoked `-latestTime` by driverFOAM's own tet workflow DAG) with `snGrad corrected` (the case's committed default).

Replaces the former `setup/mesh/tet/run_parallel_interface_sweep.sh` bash script — mesh generation, `checkMesh`, and the tet electroProperties/fvSchemes overlay swap are now handled directly by driverFOAM's own `manufacturedBathBidomain` tet workflow DAG, not by hand-rolled bash.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec setup/studies/interfaceCurrentConvergence/sweep_tet_unweightedHarmonic.json
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec setup/studies/interfaceCurrentConvergence/sweep_tet_distanceWeightedHarmonic.json
```

Each method is its own spec (rather than a single sweep with a `method` axis) because the sweep engine's `caseId`/`output_dir_name` templating can only reference scalar/list independent axes; `electro_property_overrides` (the mechanism that sets `interfaceConductivityInterpolation`) is a per-case dict and isn't safe to reference there, so it's fixed in each spec's `base` instead. Same pattern this tutorial already uses for `sweep_hex_electrodePair.json`/`sweep_hex_groundElectrode.json`.

## Status

Verified with a real `driverFoam sweep-run` at `N=10` (`unweightedHarmonic`, 2026-08-19): the case meshes, solves, and writes `bathBidomainInterfaceMetrics.csv` under the sweep's own archive layout (`<case_root>/<caseId>/setup/studies/interfaceCurrentConvergence/results/sweepCases...`, see the coupling study's README for how that layout actually works). `N=20,40,80` and the `distanceWeightedHarmonic` spec are unrun but use the identical mechanism — no reason to expect them to behave differently. Not yet checked against `reference/bath_tet_convergence.csv` numerically. `bath_predictor_corrector: true` in `base` matches `setup/mesh/tet/electroProperties`'s own baked-in default (the bash script never touched that key, so it always ran with predictor-corrector coupling enabled).

This also required a fix: the checked-in `constant/electroProperties` and `setup/mesh/tet/electroProperties` were both missing the `bidomainSolverCoeffs.{verificationModel,manufacturedBidomain}.fdaBathVariant` key that `_apply_case` always writes — every driverFOAM sweep for this tutorial (tet or hex, old specs included) crashed with a `KeyError` before this was added.

## Tracking & Outputs

All generated outputs land under the local `results/` folder, gitignored. Do not commit generated OpenFOAM data.
