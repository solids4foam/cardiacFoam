# tetTemporalControl - bidomain

## Purpose

The primary generic-Delaunay `sbdf2` ladders use `dt ~ h^2`. This four-case
control holds mesh, field equations, gradient reconstruction, linear and ODE
controls, and end time fixed while halving `dt` at `N=40` and `N=80` for both
gradient reconstructions. It tests whether the observed tet field separation is
spatially dominated over the accepted levels.

## Execution

```bash
driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/tetTemporalControl/sweep_tet_dt_half.json \
    --output-dir .tmp/driverfoam/bidomain-tet-dt-half
```
