# odeToleranceControl - bidomain

## Purpose

This two-case control repeats the finest fixed-grid temporal case in 1D and
2D with tighter RKF45 tolerances (`absTol=1e-12`, `relTol=1e-10`). Compare its
field metrics with the matching baseline temporal result, whose explicit
controls are `1e-10` and `1e-8`. It bounds adaptive-ODE sensitivity only; it
does not replace the fixed-grid temporal ladder.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/odeToleranceControl/sweep_ode_tolerance.json \
    --output-dir .tmp/driverfoam/bidomain-ode-control
```
