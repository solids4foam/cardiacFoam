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

**Corrected 2026-09-26 (tutorials-are-pointers, 5.4b-B):** as `corrector/`'s
own README now says, the two `nonorth1 (0)`/`combined (1)`-carrying cases
per pair preview cleanly against the `manufacturedBidomain` record but
cannot be committed by a real `sweep-run` yet: `PIMPLE.nNonOrthogonalCorrectors`
is absent from this case's `system/fvSolution`, and a tutorial record's
direct study-key channel writes with `add_if_missing=False`.
