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
driverFoam sweep-run --spec setup/studies/gradientScheme/sweep_current.json
driverFoam sweep-run --spec setup/studies/gradientScheme/sweep_gaussLinear.json
driverFoam sweep-run --spec setup/studies/gradientScheme/sweep_limitedCorrection.json
driverFoam sweep-run --spec setup/studies/gradientScheme/sweep_orthogonalControl.json
```

`driverFoam` is the external orchestration add-on (not part of this repo;
see the root `CLAUDE.md`).

## Status

`current` and `limitedCorrection` (`N=20`) run to completion; the
resulting `system/fvSchemes` carries the intended `default leastSquares` /
`Gauss linear limited 0.5` / `limited 0.5` triple for `limitedCorrection`,
confirming `fv_scheme_overrides` actually lands. `gaussLinear` and
`orthogonalControl` use the identical mechanism and have not been run.
See `setup/studies/coupling/README.md` for where sweep-run output actually
lands if you go on to aggregate these.

`constant/electroProperties` must set
`bidomainSolverCoeffs.{verificationModel,manufacturedBidomain}.fdaBathVariant`
— `_apply_case` always writes this key, and a driverFOAM sweep for this
tutorial (tet or hex) fails with `KeyError` without it.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, gitignored. Do not commit generated OpenFOAM data.
