# odeToleranceControl - monodomainPseudoECG

## Purpose

This three-case control bounds the adaptive RKF45 contribution at the finest
fixed-mesh temporal level in 1D, 2D, and 3D.  It repeats the temporal setup
with `absTol=1e-12` and `relTol=1e-10`, compared with the documented baseline
of `1e-10` and `1e-8`.  It is a sensitivity check, not another convergence
ladder.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/odeToleranceControl/sweep_ode_tolerance.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-ode-control
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/odeToleranceControl/sweep_ode_tolerance.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-ode-control
```

Compare each result with its matching finest temporal case in
`temporalConvergence`; retain the field and every electrode's pseudo-ECG
metric in the comparison.
