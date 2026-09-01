# linearToleranceControl - bidomain

## Purpose

The primary tetrahedral ladder solves both `Vm` and the gauge-constrained
`phiE` block with absolute tolerance `1e-15`. This four-case sensitivity
control changes only the `phiE|phiI` block to `1e-6` at the two finest levels
of both gradient reconstructions; the `Vm` block remains tight. Compare
volume-weighted and maximum norms against the matching primary cases; residual
tolerance is a monitor, not an error estimate.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bidomain/setup/studies/linearToleranceControl/sweep_tet_phi_tolerance.json \
    --output-dir .tmp/driverfoam/bidomain-linear-tolerance
```
