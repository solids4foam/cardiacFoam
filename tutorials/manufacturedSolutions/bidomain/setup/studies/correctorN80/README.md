# correctorN80 - bidomain

## Purpose

This study runs the four `(nOuterCorrectors, nNonOrthogonalCorrectors)`
combinations at `N=80`, on the generic-Delaunay mesh with tight linear and
ODE controls and a 144-step window. It is an iterative-coupling sensitivity control, testing solver-loop
convergence rather than spatial resolution.

## Execution

```bash
driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/correctorN80/sweep_corrector_n80.json \
    --output-dir .tmp/driverfoam/bidomain-corrector-n80
```

**Corrected 2026-09-26, then resolved the same day (controller decision), as
`corrector/`'s own README now says:** `system/fvSolution` now states
`nNonOrthogonalCorrectors 0;` explicitly (the default it was silently
relying on before), proven behaviour-neutral by
`regression/regressionTest.sh`. All four cases here now `strict_plan`
cleanly.
