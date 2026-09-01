# tetTemporalControl - bathBidomain

## Purpose

The distance-weighted-harmonic tetrahedral interface-current sequence uses
`dt ~ h^2` with least-squares reconstruction and predictor--corrector bath
coupling. They retain the Paper I `groundElectrode` boundary variant. These two
controls hold the mesh, interface treatment, linear tolerance, end time, and
coupling fixed while halving `dt` at `N=40` and `N=80`. They test whether
temporal/splitting error contributes to the solved interface-current behaviour,
including the finest-level anomaly.

This is a temporal-sensitivity control, not an unstructured temporal-order
claim: the dedicated Cartesian fixed-mesh ladder supplies that measurement.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/bathBidomain/setup/studies/tetTemporalControl/sweep_tet_dt_half.json \
    --output-dir .tmp/driverfoam/bathBidomain-tet-dt-half
```
