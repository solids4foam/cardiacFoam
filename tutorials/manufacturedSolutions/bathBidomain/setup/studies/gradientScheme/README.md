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

Each variant is its own sweep spec, with its schemes fixed in `base` as
`system/fvSchemes:` keys. A spec states only the entries that differ from
the case's own `system/fvSchemes`, which is `current`, so
`sweep_current.json` states none.

## Execution

Specs: `setup/studies/gradientScheme/sweep_current.json`, `setup/studies/gradientScheme/sweep_gaussLinear.json`, `setup/studies/gradientScheme/sweep_limitedCorrection.json`, `setup/studies/gradientScheme/sweep_orthogonalControl.json`.

```bash
[omnidriver command to run]
```


## Status

`current` and `limitedCorrection` (`N=20`) have been run and land the intended
scheme keys; `gaussLinear` and `orthogonalControl` use the same mechanism and
have not been run.
See `setup/studies/coupling/README.md` for where sweep-run output actually
lands if you go on to aggregate these.

These specs name no boundary variant, so they run the case's own
`electrodePair`.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, gitignored. Do not commit generated OpenFOAM data.
