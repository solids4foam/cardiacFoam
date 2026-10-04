# correctorN80 - bidomain

## Purpose

This study runs the four `(nOuterCorrectors, nNonOrthogonalCorrectors)`
combinations at `N=80`, on the generic-Delaunay mesh with tight linear and
ODE controls and a 144-step window. It is an iterative-coupling sensitivity control, testing solver-loop
convergence rather than spatial resolution.

## Execution

Spec: `tutorials/manufacturedSolutions/bidomain/setup/studies/correctorN80/sweep_corrector_n80.json`.

```bash
[omnidriver command to run]
```

`system/fvSolution` states `nNonOrthogonalCorrectors 0;` explicitly, as
`corrector/`'s README explains.
