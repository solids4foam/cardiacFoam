# temporalConvergence - bathBidomain

## Purpose

This eight-case fixed-mesh time-step ladder measures the complete configured
bath-bidomain advance, including its predictor--corrector interface workflow.
It uses the explicitly archived `godunov` coupling at `N=640` in 1D and 2D,
with four successive halvings of `dt` and fixed end time 0.02. It does not
transfer temporal order from standalone monodomain or bidomain runs.

The companion tetrahedral `tetTemporalControl` repeats selected interface
levels at `dt/2`; it asks whether temporal error contaminates the interface
current sequence rather than establishing a separate temporal order on the
unstructured mesh.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bathBidomain/setup/studies/temporalConvergence/sweep_hex_temporal_godunov.json \
    --output-dir .tmp/driverfoam/bathBidomain-temporal
```

## Status

The manifest is driverFOAM-plan validated. Numerical results and any temporal
order claim remain pending an OpenFOAM run and a check for a spatial-error
floor at the accepted levels.
