# corrector - bidomain

## Purpose

This study validates the predictor-corrector inner loop convergence and stability for the bidomain solver.

## Execution

```bash
[omnidriver command to run]
```

Sweeps `N = 10, 20, 40` across the four reported variants (`baseline`,
`outer2`, `nonorth1`, `combined`), each a `(nOuterCorrectors,
nNonOrthogonalCorrectors)` pair applied as direct
`system/fvSolution:PIMPLE.nOuterCorrectors`/`PIMPLE.nNonOrthogonalCorrectors`
study keys against the `manufacturedBidomain` tutorial record; the short,
fixed step-count screening window (2/9/36 steps) is set via the direct keys
`system/controlDict:writeControl`/`writeInterval`/`writeFormat`.

Comparing the 12 results across variants/resolutions is handled by the
postprocessing module.

Each case's raw output lives under the sweep's own `output_dir`.

`system/fvSolution` states `nNonOrthogonalCorrectors 0;` explicitly: a
study's direct keys change keys that exist and never invent one, so the
`nonorth1`/`combined` variants need it present.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
