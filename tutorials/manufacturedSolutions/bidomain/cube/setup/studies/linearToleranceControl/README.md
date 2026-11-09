# linearToleranceControl - bidomain

The [tetConvergence](../tetConvergence/README.md) ladder solves `phiE` and
`Vm` to the absolute tolerance `1e-15`. This control repeats the `N = 40`
tetrahedral level of both gradient reconstructions with the `phiE` tolerance
loosened to `1e-6` (`-phiTolerance 1e-6`), the `Vm` tolerance unchanged, and
compares the volume-weighted and maximum norms. A residual tolerance is a
stopping rule, not an error estimate; the field errors show whether the
algebraic error stays below the discretisation error.

## Ladder

| N | phiE tolerance |
|---|---|
| 40 | 1e-6, 1e-15 |

for `leastSquares` and `Gauss linear`, with the other settings of
tetConvergence, including the `Allrun` default time step
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which ends exactly at
`endTime 0.2`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_tol<tolerance>_tet_N40`, the names that `../summarize.py` reads;
the `1e-15` runs are links to the tetConvergence runs, which use that
tolerance. From this folder, after
[tetConvergence](../tetConvergence/README.md):

```bash
mkdir -p results
for g in leastSquares GaussLinear
do
    ln -s ../../tetConvergence/results/${g}_tet_N40 results/${g}_tol1e-15_tet_N40
    r=results/${g}_tol1e-6_tet_N40
    rsync -a --exclude setup ../../../ $r/
    $r/Allrun parallel -mesh tet -gradient $g -N 40 -phiTolerance 1e-6
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_phi_tolerance.json`.

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`; gauge-shifted potentials, largest `phiE` iteration count
and final residual over the run.

| gradient | phiE tolerance | Vm L2 | Vm Linf | phiE L2 | phiE Linf | phiI L2 | phiE iterations | phiE final residual |
|---|---|---|---|---|---|---|---|---|
| `leastSquares` | 1e-15 | 6.3151e-04 | 5.7885e-03 | 4.4747e-04 | 4.1461e-03 | 1.8553e-04 | 494 | 1.00e-15 |
| `leastSquares` | 1e-6 | 6.3147e-04 | 5.7891e-03 | 4.4760e-04 | 4.1435e-03 | 1.8586e-04 | 233 | 1.00e-06 |
| `leastSquares` | change | -0.01 % | +0.01 % | +0.03 % | -0.06 % | +0.18 % |  |  |
| `Gauss linear` | 1e-15 | 5.0344e-03 | 4.5213e-02 | 3.5670e-03 | 3.1940e-02 | 1.4680e-03 | 493 | 1.00e-15 |
| `Gauss linear` | 1e-6 | 5.0345e-03 | 4.5213e-02 | 3.5673e-03 | 3.1938e-02 | 1.4679e-03 | 234 | 1.00e-06 |
| `Gauss linear` | change | +0.00 % | +0.00 % | +0.01 % | -0.01 % | -0.01 % |  |  |

At `N = 40` the loose tolerance halves the `phiE` iterations and changes the
errors by at most 0.18 % (`phiI` L2, `leastSquares`); with `Gauss linear` by
at most 0.01 %. The algebraic error at `1e-6` stays below the discretisation
error at this level. Every solve reaches its tolerance below the explicit
`maxIter 5000`.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
