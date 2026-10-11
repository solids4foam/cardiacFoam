# tetTemporalControl - bidomain

The [tetConvergence](../tetConvergence/README.md) ladder refines the time
step with the mesh (`dt ~ h^2`) under first-order `godunov` coupling, so its
slopes are space-time slopes. This control keeps the mesh and halves the time
step at the two finest levels; a change much smaller than the change between
successive levels shows that the tetrahedral errors are spatial.

## Ladder

| N | dt | dt/2 |
|---|---|---|
| 40 | 0.000560538 | 0.000280269 |
| 80 | 0.000140174 | 7.0087e-05 |

for both gradient reconstructions, with the other settings of
tetConvergence.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_tet_N<n>_dt<dt>`, the names that `../summarize.py` reads. From
this folder:

```bash
mkdir -p results
for spec in "40 0.000560538 0.000280269" "80 0.000140174 7.0087e-05"
do
    set -- $spec; n=$1
    for dt in $2 $3
    do
        for g in leastSquares GaussLinear
        do
            r=results/${g}_tet_N${n}_dt$dt
            rsync -a --exclude setup ../../../ $r/
            $r/Allrun parallel -mesh tet -gradient $g -N $n -dt $dt
        done
    done
done
python3 ../summarize.py results
```

## Measured behaviour

No results yet.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
