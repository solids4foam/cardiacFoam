# tetTemporalControl - bathBidomain

The [coupling](../coupling/README.md) and
[interfaceCurrentConvergence](../interfaceCurrentConvergence/README.md)
ladders refine the time step with the mesh (`dt ~ h^2`). This control keeps
the tetrahedral mesh, the predictor-corrector, the `distanceWeightedHarmonic`
interface conductivity and `endTime 0.02`, and halves the time step at the
two finest levels; a small change shows that the interface diagnostics there
are spatial.

## Ladder

| N | dt | dt/2 |
|---|---|---|
| 40 | 0.000560538 | 0.000280269 |
| 80 | 0.0001401345 | 7.006725e-05 |

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`tet_N<n>_dt<dt>`, the names that `../summarize.py` reads (each `tet_N<n>` is
one temporal pair). From this folder:

```bash
mkdir -p results
for spec in "40 0.000560538 0.000280269" "80 0.0001401345 7.006725e-05"
do
    set -- $spec; n=$1
    for dt in $2 $3
    do
        r=results/tet_N${n}_dt$dt
        rsync -a --exclude setup ../../../ $r/
        $r/Allrun parallel -mesh tet -N $n -dt $dt -endTime 0.02
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_dt_half.json`.

## Measured behaviour

This study has not been run (`results/` is empty).

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
