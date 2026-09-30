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

Each variant is its own sweep spec: `fv_scheme_overrides` is a per-case
list-of-dicts, which the sweep engine's case-id templating cannot reference,
so it is fixed in each spec's `base` rather than swept as an axis.

## Execution

Specs: `setup/studies/gradientScheme/sweep_current.json`, `setup/studies/gradientScheme/sweep_gaussLinear.json`, `setup/studies/gradientScheme/sweep_limitedCorrection.json`, `setup/studies/gradientScheme/sweep_orthogonalControl.json`.

```bash
[omnidriver command to run]
```

`omnidriver` is the external orchestration add-on (not part of this repo;
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
— `_apply_case` always writes this key, and an omnidriver sweep for this
tutorial (tet or hex) fails with `KeyError` without it.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, gitignored. Do not commit generated OpenFOAM data.
