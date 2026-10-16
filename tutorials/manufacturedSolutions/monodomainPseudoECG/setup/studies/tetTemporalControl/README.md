# tetTemporalControl - monodomainPseudoECG

The [tetConvergence](../tetConvergence/README.md) ladder refines the time
step with the mesh (`dt ~ h^2`), so its slopes are space-time slopes. This
control keeps the mesh and halves the time step at the two finest levels; a
small change shows that the tetrahedral errors are spatial.

## Ladder

| N | dt | dt/2 |
|---|---|---|
| 40 | 0.000560538 | 0.000280269 |
| 80 | 0.000140174 | 7.0087e-05 |

for both gradient reconstructions and both fields, as in tetConvergence.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_<field>_tet_N<n>_dt<dt>`, the names that `../summarize.py` reads.
From this folder:

```bash
mkdir -p results
for spec in "40 0.000560538 0.000280269" "80 0.000140174 7.0087e-05"
do
    set -- $spec; n=$1
    for dt in $2 $3
    do
        for g in leastSquares GaussLinear
        do
            for f in fda anisotropic
            do
                r=results/${g}_${f}_tet_N${n}_dt$dt
                rsync -a --exclude setup ../../../ $r/
                $r/Allrun parallel -mesh tet -gradient $g -field $f -N $n -dt $dt
            done
        done
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_dt_half.json`.

## Measured behaviour

Results of the last run of this study (`results/`), `leastSquares`,
anisotropic field, `N = 40`:

| field | dt | L2 | Linf |
|---|---|---|---|
| Vm | 0.000560538 | 1.74622e-04 | 1.65335e-03 |
| Vm | 0.000280269 | 1.74580e-04 | 1.64906e-03 |
| u1 | 0.000560538 | 4.89273e-05 | 4.34089e-04 |
| u1 | 0.000280269 | 4.78638e-05 | 4.50768e-04 |

Halving the time step changes `Vm` by 0.024 % in L2 and 0.26 % in Linf, two
orders below the change between successive tetrahedral levels: the
tetrahedral `Vm` errors are spatial.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
