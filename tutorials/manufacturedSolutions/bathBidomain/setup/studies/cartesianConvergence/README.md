# cartesianConvergence - bathBidomain

Spatial convergence on Cartesian hexahedra (`Allrun`, `-mesh hex`) in 1D, 2D
and 3D: the ground-electrode problem with the predictor-corrector bath
coupling (the case default) and without it (`-predictorCorrector no`), the
electrode-pair variant, and the assembled interface current of the exact
`phiE` field on the 3D boxes, where every face is orthogonal and the
conductivity diagonal.

## Ladder

| N | dt | hex cells (3D) |
|---|---|---|
| 10 | 0.00892857 | 3000 |
| 20 | 0.00224215 | 24000 |
| 40 | 0.000560538 | 192000 |
| 80 | 0.000140174 | 1536000 |

The time step falls as h^2, so the rates are combined space-time rates.
These steps do not divide `endTime` (0.2, or 0.02 for the electrode pair and
the exact-field control), so each run stops at the last step before it (the
verifier writes the final time in its summary and evaluates the exact field
there). Without `-dt` the case uses
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which ends exactly at
`endTime`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<configuration>_<d>D_N<n>`, the names that `../summarize.py` reads. From
this folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.00892857 0.00224215 0.000560538 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    for d in 1 2 3
    do
        rsync -a --exclude setup ../../../ results/predictor_${d}D_N$n/
        results/predictor_${d}D_N$n/Allrun parallel -dim $d -N $n -dt $dt
        rsync -a --exclude setup ../../../ results/baseline_${d}D_N$n/
        results/baseline_${d}D_N$n/Allrun parallel -dim $d -N $n -dt $dt -predictorCorrector no
        rsync -a --exclude setup ../../../ results/electrodePair_${d}D_N$n/
        results/electrodePair_${d}D_N$n/Allrun parallel -dim $d -N $n -dt $dt -endTime 0.02 -variant electrodePair
    done
    rsync -a --exclude setup ../../../ results/exactInterface_3D_N$n/
    results/exactInterface_3D_N$n/Allrun parallel -N $n -dt $dt -endTime 0.02
done
python3 ../summarize.py results
```

The omnidriver specs for this study are `sweep_hex_baseline.json`, `sweep_hex_convergence.json`, `sweep_hex_electrodePair.json` and `sweep_hex_exactInterface.json`.

The exact-field rows are in `bathBidomainInterfaceMetricsExactField.csv` of
every run; the `exactInterface` runs evaluate them at `t = 0.02`.

## Measured behaviour

Results of the last run of this study (`results/`).

Ground electrode with the predictor-corrector (`predictor`, the case
default), L2 errors:

| dim | N | Vm | order | phiE | order | phiI | order |
|---|---|---|---|---|---|---|---|
| 1D | 10 | 2.7785e-04 | | 1.2531e-02 | | 9.4920e-03 | |
| 1D | 20 | 6.9473e-05 | 2.00 | 3.1383e-03 | 2.00 | 2.3749e-03 | 2.00 |
| 1D | 40 | 1.7327e-05 | 2.00 | 7.8432e-04 | 2.00 | 5.9329e-04 | 2.00 |
| 1D | 80 | 4.3262e-06 | 2.00 | 1.9557e-04 | 2.00 | 1.4783e-04 | 2.00 |
| 2D | 10 | 3.2243e-04 | | 1.2570e-02 | | 9.4734e-03 | |
| 2D | 20 | 8.1030e-05 | 1.99 | 3.1492e-03 | 2.00 | 2.3709e-03 | 2.00 |
| 2D | 40 | 2.0235e-05 | 2.00 | 7.8712e-04 | 2.00 | 5.9236e-04 | 2.00 |
| 2D | 80 | 5.0541e-06 | 2.00 | 1.9628e-04 | 2.00 | 1.4760e-04 | 2.00 |
| 3D | 10 | 3.4466e-04 | | 1.2591e-02 | | 9.4652e-03 | |
| 3D | 20 | 8.6781e-05 | 1.99 | 3.1550e-03 | 2.00 | 2.3691e-03 | 2.00 |
| 3D | 40 | 2.1683e-05 | 2.00 | 7.8862e-04 | 2.00 | 5.9194e-04 | 2.00 |
| 3D | 80 | 5.4161e-06 | 2.00 | 1.9666e-04 | 2.00 | 1.4750e-04 | 2.00 |

Interface rows of the same runs in 3D (L2; the current-density jump and the
intracellular leakage coincide on these meshes):

| N | heart phiE | bath phiE | x=0 current jump | x=1 current jump | order (jump) |
|---|---|---|---|---|---|
| 10 | 9.6553e-03 | 1.3827e-02 | 3.3875e-03 | 3.3704e-03 | |
| 20 | 2.4176e-03 | 3.4652e-03 | 4.3265e-04 | 4.3201e-04 | 2.97 |
| 40 | 6.0413e-04 | 8.6622e-04 | 5.4361e-05 | 5.4335e-05 | 2.99 |
| 80 | 1.5045e-04 | 2.1591e-04 | 6.8023e-06 | 6.8011e-06 | 3.00 |

The heart and bath `phiE` rows converge at second order in every dimension.
The assembled current density of the solved field stays at the round-off
level of the conservative face flux (1e-9 to 1e-6) and has no order.

Ground electrode without the predictor-corrector (`baseline`), L2 errors:

| dim | N | Vm | order | phiE | order | phiI | order |
|---|---|---|---|---|---|---|---|
| 1D | 10 | 3.4424e-04 | | 1.6353e-02 | | 1.2566e-02 | |
| 1D | 20 | 8.6496e-05 | 1.99 | 4.1091e-03 | 1.99 | 3.1547e-03 | 1.99 |
| 1D | 40 | 2.1601e-05 | 2.00 | 1.0284e-03 | 2.00 | 7.8931e-04 | 2.00 |
| 1D | 80 | 5.3974e-06 | 2.00 | 2.5743e-04 | 2.00 | 1.9750e-04 | 2.00 |
| 2D | 10 | 3.8657e-04 | | 1.6391e-02 | | 1.2546e-02 | |
| 2D | 20 | 9.7501e-05 | 1.99 | 4.1201e-03 | 1.99 | 3.1503e-03 | 1.99 |
| 2D | 40 | 2.4373e-05 | 2.00 | 1.0313e-03 | 2.00 | 7.8830e-04 | 2.00 |
| 2D | 80 | 6.0910e-06 | 2.00 | 2.5816e-04 | 2.00 | 1.9725e-04 | 2.00 |
| 3D | 10 | 4.0811e-04 | | 1.6412e-02 | | 1.2536e-02 | |
| 3D | 20 | 1.0310e-04 | 1.98 | 4.1260e-03 | 1.99 | 3.1483e-03 | 1.99 |
| 3D | 40 | 2.5782e-05 | 2.00 | 1.0328e-03 | 2.00 | 7.8783e-04 | 2.00 |
| 3D | 80 | 6.4435e-06 | 2.00 | 2.5854e-04 | 2.00 | 1.9714e-04 | 2.00 |

The predictor-corrector lowers the `phiE` and `phiI` errors by 23 to 25 %
and the `Vm` error by 16 to 20 % at every level; both are second order.

Electrode pair (`-endTime 0.02`), L2 errors after removing the volume means
of `phiE` and `phiI`:

| dim | N | final time | Vm | order | phiE | order | phiI | order |
|---|---|---|---|---|---|---|---|---|
| 1D | 10 | 0.0178571 | 3.4857e-05 | | 7.1966e-03 | | 7.1977e-06 | |
| 1D | 20 | 0.0201793 | 8.2772e-06 | 2.07 | 1.8030e-03 | 2.00 | 1.4638e-06 | 2.30 |
| 1D | 40 | 0.0201794 | 1.9728e-06 | 2.07 | 4.5108e-04 | 2.00 | 3.0024e-07 | 2.29 |
| 1D | 80 | 0.0200449 | 4.8438e-07 | 2.03 | 1.1304e-04 | 2.00 | 5.4850e-07 | -0.87 |
| 2D | 10 | 0.0178571 | 3.8936e-05 | | 7.1993e-03 | | 8.1182e-06 | |
| 2D | 20 | 0.0201793 | 9.5217e-06 | 2.03 | 1.8039e-03 | 2.00 | 1.7437e-06 | 2.22 |
| 2D | 40 | 0.0201794 | 2.2893e-06 | 2.06 | 4.5130e-04 | 2.00 | 3.1085e-07 | 2.49 |
| 2D | 80 | 0.0200449 | 5.6332e-07 | 2.02 | 1.1310e-04 | 2.00 | 5.3776e-07 | -0.79 |
| 3D | 10 | 0.0178571 | 4.1010e-05 | | 7.2008e-03 | | 8.5598e-06 | |
| 3D | 20 | 0.0201793 | 1.0143e-05 | 2.02 | 1.8043e-03 | 2.00 | 1.8712e-06 | 2.19 |
| 3D | 40 | 0.0201794 | 2.4465e-06 | 2.05 | 4.5142e-04 | 2.00 | 3.1340e-07 | 2.58 |
| 3D | 80 | 0.0200449 | 6.0248e-07 | 2.02 | 1.1313e-04 | 2.00 | 5.3167e-07 | -0.76 |

`Vm` and `phiE` are second order; `phiI`, three to four orders below
`phiE`, stops falling at about 5e-7 at `N = 80`.

Assembled interface current density of the exact `phiE` field, 3D,
`-endTime 0.02`:

| N | final time | x=0 | x=1 | order | x=0 current jump | order |
|---|---|---|---|---|---|---|
| 10 | 0.0178571 | 2.7032e-03 | 2.7032e-03 | | 3.1169e-03 | |
| 20 | 0.0201793 | 1.3552e-03 | 1.3552e-03 | 1.00 | 3.9860e-04 | 2.97 |
| 40 | 0.0201794 | 6.7787e-04 | 6.7787e-04 | 1.00 | 5.0108e-05 | 2.99 |
| 80 | 0.0200449 | 3.3895e-04 | 3.3895e-04 | 1.00 | 6.2846e-06 | 3.00 |

On these orthogonal faces the two-point flux of the exact field is first
order at both interfaces: the curvature of the exact field leaves an `O(h)`
error between the cell-centre difference and the face current.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
