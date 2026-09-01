# correctorN80 - bidomain

## Purpose

This endpoint companion to `corrector` repeats the four
`(nOuterCorrectors, nNonOrthogonalCorrectors)` combinations at `N=80`. It
closes the fine-level gap left by the existing `N=10,20,40` screen, using the
same generic-Delaunay mesh, tight linear and ODE controls, and a 144-step
window. It is an iterative-coupling sensitivity control, not an additional
spatial ladder.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/correctorN80/sweep_corrector_n80.json \
    --output-dir .tmp/driverfoam/bidomain-corrector-n80
```
