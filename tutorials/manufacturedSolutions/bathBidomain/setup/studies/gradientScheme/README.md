# gradientScheme - bathBidomain

## Purpose

Screens four gradient/laplacian/snGrad scheme combinations at `N=20` on the
conformal tetrahedral mesh (`distanceWeightedHarmonic` interface method,
fixed):

| variant | gradSchemes.default | laplacianSchemes.default | snGradSchemes.default |
|---|---|---|---|
| `current` (checked-in default) | `leastSquares` | `Gauss linear corrected` | `corrected` |
| `gaussLinear` | `Gauss linear` | `Gauss linear corrected` | `corrected` |
| `limitedCorrection` | `leastSquares` | `Gauss linear limited 0.5` | `limited 0.5` |
| `orthogonalControl` | `leastSquares` | `Gauss linear orthogonal` | `orthogonal` |

Formerly `run_gradient_scheme_screen.sh`, which delegated to
`run_parallel_interface_sweep.sh` per variant. Each variant is now its own
sweep spec (`fv_scheme_overrides` is a per-case list-of-dicts, which the
sweep engine's case-id templating can't safely reference, so it's fixed in
each spec's `base` rather than swept as an axis) — same reasoning as
`setup/studies/interfaceCurrentConvergence/`.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec setup/studies/gradientScheme/sweep_current.json
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec setup/studies/gradientScheme/sweep_gaussLinear.json
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec setup/studies/gradientScheme/sweep_limitedCorrection.json
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec setup/studies/gradientScheme/sweep_orthogonalControl.json
```

## Status

Verified with real `driverFoam sweep-run`s at `N=20` for `current` and
`limitedCorrection` (2026-08-19): both complete, and the resulting
`system/fvSchemes` carries the intended `default leastSquares` /
`Gauss linear limited 0.5` / `limited 0.5` triple for `limitedCorrection`,
confirming `fv_scheme_overrides` actually lands. `gaussLinear` and
`orthogonalControl` unrun in this session but use the identical mechanism.
See `setup/studies/coupling/README.md` for where sweep-run output actually
lands (not the naive `archive_dir_name` reading) if you go on to aggregate
these.

This also required a fix: the checked-in `constant/electroProperties` and
`setup/mesh/tet/electroProperties` were both missing the
`bidomainSolverCoeffs.{verificationModel,manufacturedBidomain}.fdaBathVariant`
key that `_apply_case` always writes — every driverFOAM sweep for this
tutorial (tet or hex, old specs included) crashed with a `KeyError` before
this was added.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, gitignored. Do not commit generated OpenFOAM data.
