# gradientScheme - bathBidomain

Cell-gradient reconstruction and non-orthogonal correction on the conformal
tetrahedral `N = 20` mesh of [tetConvergence](../tetConvergence/README.md),
with the predictor-corrector and `distanceWeightedHarmonic` interface
conductivity. Each variant changes the defaults of `gradSchemes`,
`laplacianSchemes` and `snGradSchemes`, for `Vm` and `phiE` alike:

| variant | `Allrun` options | gradSchemes | laplacianSchemes | snGradSchemes |
|---|---|---|---|---|
| `current` | (defaults) | `leastSquares` | `Gauss linear corrected` | `corrected` |
| `gaussLinear` | `-gradient GaussLinear` | `Gauss linear` | `Gauss linear corrected` | `corrected` |
| `limited` | `-snGrad limited` | `leastSquares` | `Gauss linear limited 0.5` | `limited 0.5` |
| `orthogonal` | `-snGrad orthogonal` | `leastSquares` | `Gauss linear orthogonal` | `orthogonal` |

## Ladder

One level: `N = 20` (111482 cells with gmsh 4.15.2, seed 1),
`dt = 0.00224215`, `endTime 0.02`; the step does not divide `endTime`, so
each run stops at 0.0201793.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<variant>_tet_N20`, the names that `../summarize.py` reads. From this
folder:

```bash
mkdir -p results
for spec in "current" "gaussLinear -gradient GaussLinear" \
            "limited -snGrad limited" "orthogonal -snGrad orthogonal"
do
    set -- $spec; v=$1; shift
    rsync -a --exclude setup ../../../ results/${v}_tet_N20/
    results/${v}_tet_N20/Allrun parallel -mesh tet -N 20 -dt 0.00224215 -endTime 0.02 "$@"
done
python3 ../summarize.py results
python3 ../summarize.py results --relative-to current_tet
```

The omnidriver specs for this study are `sweep_current.json`, `sweep_gaussLinear.json`, `sweep_limitedCorrection.json` and `sweep_orthogonalControl.json`.

## Measured behaviour

Results of the last run of this study (`results/`), L2:

| variant | Vm | phiE | phiI | heart phiE | bath phiE | x=0 current jump | x=0 leakage | x=0 assembled | x=1 assembled |
|---|---|---|---|---|---|---|---|---|---|
| current | 1.5160e-04 | 1.2795e-03 | 9.8885e-04 | 9.9786e-04 | 1.3992e-03 | 5.1075e-04 | 6.6865e-04 | 2.6704e-04 | 2.8760e-04 |
| gaussLinear | 1.5673e-03 | 3.2652e-02 | 2.9378e-02 | 2.9405e-02 | 3.4160e-02 | 1.2565e-03 | 7.8893e-04 | 8.0624e-04 | 8.2557e-04 |
| limited | 1.1883e-03 | 8.6495e-04 | 5.1036e-04 | 9.3496e-04 | 8.2774e-04 | 5.6849e-04 | 5.2617e-04 | 7.6972e-04 | 8.0422e-04 |
| orthogonal | 3.3298e-03 | 5.2453e-03 | 3.8394e-03 | 4.4938e-03 | 5.5833e-03 | 1.6455e-03 | 1.1496e-03 | 1.3847e-03 | 1.4042e-03 |

`Gauss linear` raises the `phiE` errors 25-fold and `Vm` 10-fold;
`orthogonal` (no non-orthogonal correction) raises every row by 1.7 to 22
times. `limited 0.5` lowers the `phiE` and `phiI` errors by 32 and 48 % but
raises `Vm` 8-fold and the assembled current 2.9-fold. The defaults give the
smallest `Vm` and assembled-current errors.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
