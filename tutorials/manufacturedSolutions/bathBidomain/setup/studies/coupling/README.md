# coupling - bathBidomain

Effect of the predictor-corrector bath coupling (`bathPredictorCorrector`,
`Allrun -predictorCorrector yes|no`) on the same conformal tetrahedral meshes
as [tetConvergence](../tetConvergence/README.md). Without it each step solves
the global `phiE` and then `Vm` with that `phiE` held fixed; with it each step
predicts `Vm`, solves `phiE` and corrects `Vm`. This is a same-mesh
sensitivity, not a convergence study.

## Ladder

| N | dt | endTime |
|---|---|---|
| 10 | 0.00892857 | 0.02 |
| 20 | 0.00224215 | 0.02 |
| 40 | 0.000560538 | 0.02 |
| 80 | 0.000140174 | 0.02 |

for `-predictorCorrector yes` (`predictor`) and `no` (`baseline`). The steps
do not divide `endTime`; each run stops at the last step before it.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<configuration>_tet_N<n>`, the names that `../summarize.py` reads. From this
folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.00892857 0.00224215 0.000560538 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    for pc in yes no
    do
        [ $pc = yes ] && c=predictor || c=baseline
        rsync -a --exclude setup ../../../ results/${c}_tet_N$n/
        results/${c}_tet_N$n/Allrun parallel -mesh tet -N $n -dt $dt -endTime 0.02 -predictorCorrector $pc
    done
done
python3 ../summarize.py results
python3 ../summarize.py results --relative-to baseline_tet
```

## Measured behaviour

Results of the last run of this study (`results/`, `N = 10, 20, 40`):
change of the predictor-corrector against the baseline.

| N | myocardium phiE L2 | bath phiE L2 | x=0 assembled current density L2 |
|---|---|---|---|
| 10 | -53.6 % | -50.7 % | -0.8 % |
| 20 | -46.2 % | -44.6 % | -1.3 % |
| 40 | -40.4 % | -38.9 % | +0.07 % |

Interface rows of the predictor-corrector runs (L2):

| N | heart phiE | bath phiE | x=0 current jump | x=0 leakage | x=0 assembled | x=1 assembled |
|---|---|---|---|---|---|---|
| 10 | 3.3075e-03 | 4.6624e-03 | 3.5185e-03 | 3.3821e-03 | 4.4417e-04 | 4.6773e-04 |
| 20 | 1.0816e-03 | 1.5172e-03 | 5.6127e-04 | 5.3635e-04 | 2.8138e-04 | 3.0119e-04 |
| 40 | 3.3820e-04 | 4.7245e-04 | 1.1232e-04 | 1.3520e-04 | 2.2834e-04 | 2.4016e-04 |

The predictor-corrector halves the regional `phiE` errors at `N = 10` and
lowers them by about 40 % at `N = 40`; it leaves the assembled interface
current within 1.3 %.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
