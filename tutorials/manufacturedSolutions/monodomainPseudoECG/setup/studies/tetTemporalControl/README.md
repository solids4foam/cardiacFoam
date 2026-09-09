# tetTemporalControl - monodomainPseudoECG

## Purpose

The generic-Delaunay spatial ladders use `dt ~ h^2`; they are combined
space--time paths under Lie--Godunov coupling.  This eight-case control holds
the tetrahedral mesh, conductivity tensor, gradient reconstruction, solver,
reference quadrature, and end time fixed, while halving `dt` at `N=40` and
`N=80`.  It covers all four primary configurations (axis/rotated tensor times
Gauss--linear/least-squares).

The comparison is accepted only when the fixed-mesh change is smaller than the
field-error separation used to support the reconstruction conclusion.  Archive
the pseudo-ECG metrics per electrode together with the field metrics.

## Execution

```bash
driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetTemporalControl/sweep_tet_dt_half.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-tet-dt-half
driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetTemporalControl/sweep_tet_dt_half.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-tet-dt-half
```
