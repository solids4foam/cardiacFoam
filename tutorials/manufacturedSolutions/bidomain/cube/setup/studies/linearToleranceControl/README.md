# linearToleranceControl - bidomain

The [tetConvergence](../tetConvergence/README.md) ladder solves `phiE` and
`Vm` to the absolute tolerance `1e-15`. This control repeats the two finest
tetrahedral levels of both gradient reconstructions with the `phiE` tolerance
loosened to `1e-6` (`-phiTolerance 1e-6`), the `Vm` tolerance unchanged, and
compares the volume-weighted and maximum norms. A residual tolerance is a
stopping rule, not an error estimate; the field errors show whether the
algebraic error stays below the discretisation error.

## Ladder

| N | dt | phiE tolerance |
|---|---|---|
| 40 | 0.000560538 | 1e-6, 1e-15 |
| 80 | 0.000140174 | 1e-6, 1e-15 |

for `leastSquares` and `Gauss linear`, with the other settings of
tetConvergence.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_tol<tolerance>_tet_N<n>`, the names that `../summarize.py` reads.
From this folder:

```bash
mkdir -p results
for spec in "40 0.000560538" "80 0.000140174"
do
    set -- $spec; n=$1; dt=$2
    for g in leastSquares GaussLinear
    do
        for tol in 1e-6 1e-15
        do
            r=results/${g}_tol${tol}_tet_N$n
            rsync -a --exclude setup ../../../ $r/
            $r/Allrun parallel -mesh tet -gradient $g -N $n -dt $dt -phiTolerance $tol
        done
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_phi_tolerance.json`.

## Measured behaviour

Results of the last run of this study (`results/`), gauge-shifted `phiE`;
the order is between `N = 40` and `N = 80`.

| gradient | quantity | tolerance | N = 40 | N = 80 | order |
|---|---|---|---|---|---|
| leastSquares | Vm L2 | 1e-6 | 6.3160e-04 | 1.6115e-04 | 1.97 |
| leastSquares | Vm L2 | 1e-15 | 6.4760e-04 | 1.6484e-04 | 1.97 |
| leastSquares | Vm Linf | 1e-6 | 5.7859e-03 | 1.0465e-02 | -0.85 |
| leastSquares | Vm Linf | 1e-15 | 5.7855e-03 | 1.7167e-03 | 1.75 |
| leastSquares | phiE L2 | 1e-6 | 4.4724e-04 | 1.2414e-04 | 1.85 |
| leastSquares | phiE L2 | 1e-15 | 4.5984e-04 | 1.1699e-04 | 1.97 |
| leastSquares | phiE Linf | 1e-6 | 4.1234e-03 | 1.0779e-02 | -1.39 |
| leastSquares | phiE Linf | 1e-15 | 4.1261e-03 | 1.2239e-03 | 1.75 |
| GaussLinear | Vm L2 | 1e-6 | 5.0349e-03 | 3.7599e-03 | 0.42 |
| GaussLinear | Vm Linf | 1e-6 | 4.5221e-02 | 2.7108e-02 | 0.74 |
| GaussLinear | phiE L2 | 1e-6 | 3.5643e-03 | 2.6635e-03 | 0.42 |
| GaussLinear | phiE Linf | 1e-6 | 3.1948e-02 | 1.9194e-02 | 0.74 |

With `leastSquares` the loose tolerance leaves the L2 errors at the same
scale, but at `N = 80` it raises the maximum errors of `Vm` and `phiE` six-
and nine-fold, so their orders turn negative; at `1e-15` both maximum norms
converge. The `Gauss linear` errors are dominated by the reconstruction and
the loose tolerance does not change their scale (compare
[tetConvergence](../tetConvergence/README.md)).

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
