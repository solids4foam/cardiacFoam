# weakCouplingControl - monodomain1D3D

The ladder of [coupledConvergence](../coupledConvergence/README.md) with
unidirectional explicit transfer and `R_pvj = 1e6` (`-rPvj 1e6`): the junction
current is six orders smaller, so any error the junction adds to the tissue
or the graph shows as a difference from the `R_pvj = 1` runs.

## Ladder

| N | graph | dt |
|---|---|---|
| 10 | nodes011 | 0.008971136 |
| 20 | nodes021 | 0.002242784 |
| 40 | nodes041 | 0.000560696 |
| 80 | nodes081 | 0.000140174 |

`endTime 0.1`; the final times are those of coupledConvergence.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`rPvj1e6_N<n>`, the names that `../summarize.py` reads. From this folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.008971136 0.002242784 0.000560696 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    rsync -a --exclude setup ../../../ results/rPvj1e6_N$n/
    results/rPvj1e6_N$n/Allrun -N $n -dt $dt -rPvj 1e6
done
python3 ../summarize.py results
```

## Measured behaviour

Results of the last run of this study (`results/`).

| N | h | 3D Vm L2 | 3D Vm Linf | 1D Vm L2 | 1D Vm Linf | source error L1 | order 3D L2 | order 3D Linf | order 1D L2 | order 1D Linf |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.1 | 1.8571e-03 | 4.9598e-03 | 1.8097e-04 | 3.4317e-04 | 2.3941e-10 | | | | |
| 20 | 0.05 | 4.9054e-04 | 1.3812e-03 | 4.5652e-05 | 8.8475e-05 | 1.5823e-10 | 1.92 | 1.84 | 1.99 | 1.96 |
| 40 | 0.025 | 1.2217e-04 | 3.4859e-04 | 1.1195e-05 | 2.1956e-05 | 4.3300e-11 | 2.01 | 1.99 | 2.03 | 2.01 |
| 80 | 0.0125 | 3.0648e-05 | 8.7734e-05 | 2.7873e-06 | 5.5004e-06 | 1.1445e-11 | 1.99 | 1.99 | 2.01 | 2.00 |

Against the `R_pvj = 1` unidirectional explicit runs, the 3D L2 error is
0.4% larger at every level and the 3D Linf error differs by 0.05% at `N = 10`
and not at all to the printed digits from `N = 20`. The 3D field converges at
the same order with or without the junction, so its convergence is evidence
for the monodomain operator, not for the junction; the junction is checked by
the source error of coupledConvergence. The graph error equals that of the
`R_pvj = 1` runs: under unidirectional transfer the graph receives no current
and its manufactured source has no junction term.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
