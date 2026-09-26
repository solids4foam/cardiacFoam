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
