# temporalConvergence - bathBidomain

Fixed-mesh temporal convergence of the complete bath advance (Godunov
splitting of the ionic step, predictor-corrector bath coupling, `backward`
`Vm`) on Cartesian hexahedra with `N = 640`, ground electrode, `endTime
0.02`. The mesh is fine enough that the time step dominates the `Vm` error.

## Ladder

| dimension | N | dt |
|---|---|---|
| 1D, 2D | 640 | 0.001121075, 0.000560538, 0.000280269, 0.0001401345 |

The steps do not divide `endTime 0.02`, so each run stops at the last step
before it (0.0201793, 0.0201794, 0.0198991 and 0.0200392); the verifier
evaluates the exact field at that time.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<d>D_dt<dt>`, the names that `../summarize.py` reads (each `<d>D` is one
temporal ladder). From this folder:

```bash
mkdir -p results
for d in 1 2
do
    for dt in 0.001121075 0.000560538 0.000280269 0.0001401345
    do
        rsync -a --exclude setup ../../../ results/${d}D_dt$dt/
        results/${d}D_dt$dt/Allrun parallel -dim $d -N 640 -dt $dt -endTime 0.02
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_hex_temporal_godunov.json`.

## Measured behaviour

Results of the last run of this study (`results/`, 1D only), L2:

| dt | Vm | order | phiE | order | phiI | order |
|---|---|---|---|---|---|---|
| 0.001121075 | 4.6616e-06 | | 6.6860e-06 | | 2.0117e-06 | |
| 0.000560538 | 2.2958e-06 | 1.02 | 4.3020e-06 | 0.64 | 1.6177e-06 | 0.31 |
| 0.000280269 | 1.1231e-06 | 1.03 | 3.1622e-06 | 0.44 | 1.5175e-06 | 0.09 |
| 0.0001401345 | 5.6265e-07 | 1.00 | 2.6340e-06 | 0.26 | 1.4988e-06 | 0.02 |

`Vm` is first order, the order of the Godunov splitting. `phiE` and `phiI`
reach the spatial error of the `N = 640` mesh (about 2.5e-6 and 1.5e-6)
within the ladder, so their observed orders are not temporal.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
