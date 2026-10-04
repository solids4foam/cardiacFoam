# interfaceCurrentConvergence - bathBidomain

## Purpose

Sweeps assembled interface-current rows for `@tbl-bath-bidomain-tet`: the conformal tetrahedral bath-bidomain system at `N=10,20,40,80`, `matchedSubmesh` flux assembly, comparing the `unweightedHarmonic` and `distanceWeightedHarmonic` interface-conductivity interpolation methods. Both write the assembled-current metrics computed by the `bathBidomainInterfaceMetrics` function object (invoked `-latestTime` by the omniD `manufacturedBathBidomain` record's tet route) with `snGrad corrected` (the case's committed default).

Mesh generation and `checkMesh` are the record's tet route.

## Execution

```bash
[omnidriver command to run]
```

Each method is its own spec. `sweep_tet_unweightedHarmonic.json` sets `bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation` in its `base`; `distanceWeightedHarmonic` is the case's own value, so its spec states nothing.

## Status

Run with omnidriver at `N=10` (`unweightedHarmonic`) only: the case meshes, solves, and writes `bathBidomainInterfaceMetrics.csv`. `N=20,40,80` and the `distanceWeightedHarmonic` spec use the same mechanism and are not yet run. Predictor-corrector coupling is the case's own `bathPredictorCorrector yes`.

## Tracking & Outputs

All generated outputs land under the local `results/` folder, gitignored. Do not commit generated OpenFOAM data.
