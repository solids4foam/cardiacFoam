# temporalConvergence - bathBidomain

## Purpose

This eight-case fixed-mesh time-step ladder measures the complete configured
bath-bidomain advance, including its predictor--corrector interface workflow.
It uses the explicitly archived `godunov` coupling and the Paper I
`groundElectrode` boundary variant at `N=640` in 1D and 2D, with four
successive halvings of `dt` and fixed end time 0.02. It establishes temporal
order for this configured bath-bidomain advance.

## Execution

```bash
driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bathBidomain/setup/studies/temporalConvergence/sweep_hex_temporal_godunov.json \
    --output-dir .tmp/driverfoam/bathBidomain-temporal
```

## Status

The manifest is driverFOAM-plan validated. Numerical results and any temporal
order claim remain pending an OpenFOAM run and a check for a spatial-error
floor at the accepted levels.
