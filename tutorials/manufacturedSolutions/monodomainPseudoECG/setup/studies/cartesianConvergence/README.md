# cartesianConvergence - monodomainPseudoECG

Spatial convergence on Cartesian hexahedra (`Allrun`, `-mesh hex`): the FDA
cosine field in 1D, 2D and 3D, the anisotropic sin^2 field with the full
rotated conductivity in 3D, and the 3D pseudo-ECG against its reference
integral.

## Ladder

| N | dt | hex cells (3D) |
|---|---|---|
| 10 | 0.00892857 | 1000 |
| 20 | 0.00224215 | 8000 |
| 40 | 0.000560538 | 64000 |
| 80 | 0.000140174 | 512000 |

The time step falls as h^2 so the temporal error stays below the spatial one.
These steps do not divide `endTime 0.2`, so each run stops at the last step
before 0.2 (the verifier writes the final time in its summary and evaluates
the exact field there). Without `-dt` the case uses
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which ends exactly at 0.2.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<field>_<d>D_N<n>`, the names that `../summarize.py` reads. From this
folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.00892857 0.00224215 0.000560538 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    for d in 1 2 3
    do
        rsync -a --exclude setup ../../../ results/fda_${d}D_N$n/
        results/fda_${d}D_N$n/Allrun -dim $d -field fda -N $n -dt $dt
    done
    rsync -a --exclude setup ../../../ results/anisotropic_3D_N$n/
    results/anisotropic_3D_N$n/Allrun -field anisotropic -N $n -dt $dt
done
python3 ../summarize.py results
```

The closest omnidriver spec for this study is `sweep_hex_convergence.json`: the record removes `ecgDomains` for the 1D and 2D runs, which have no ECG.

## Measured behaviour

Results of the last run of this study (`results/`).

FDA field, `Vm`:

| dim | N | L2 | Linf | order L2 | order Linf |
|---|---|---|---|---|---|
| 1D | 10 | 2.1756e-04 | 4.4091e-04 | | |
| 1D | 20 | 5.4392e-05 | 1.1114e-04 | 2.00 | 1.99 |
| 1D | 40 | 1.3559e-05 | 2.7753e-05 | 2.00 | 2.00 |
| 1D | 80 | 3.3858e-06 | 6.9334e-06 | 2.00 | 2.00 |
| 2D | 10 | 2.2887e-03 | 4.4687e-03 | | |
| 2D | 20 | 5.9466e-04 | 1.2128e-03 | 1.94 | 1.88 |
| 2D | 40 | 1.4997e-04 | 3.0912e-04 | 1.99 | 1.97 |
| 2D | 80 | 3.7526e-05 | 7.7558e-05 | 2.00 | 1.99 |
| 3D | 10 | 3.9066e-03 | 1.0271e-02 | | |
| 3D | 20 | 1.0216e-03 | 2.8541e-03 | 1.94 | 1.85 |
| 3D | 40 | 2.5804e-04 | 7.3182e-04 | 1.99 | 1.96 |
| 3D | 80 | 6.4593e-05 | 1.8387e-04 | 2.00 | 1.99 |

The auxiliary states `u1` and `u2` converge at second order in 1D and 3D; in
2D their errors are four orders below `Vm` and the orders at `N = 20`-`40`
(1.4-1.8) are pre-asymptotic, 1.9 at the finest pair.

FDA field, 3D pseudo-ECG (E1-E5, error against the order-96 reference
integral; largest and mean over the electrodes of each electrode's Linf over
time):

| N | max | order | mean | order |
|---|---|---|---|---|
| 10 | 3.7094e-04 | | 1.1939e-04 | |
| 20 | 9.4951e-05 | 1.97 | 2.9998e-05 | 1.99 |
| 40 | 2.3894e-05 | 1.99 | 7.5130e-06 | 2.00 |
| 80 | 5.9839e-06 | 2.00 | 1.8792e-06 | 2.00 |

Anisotropic field, 3D (`Vm`, and the pseudo-ECG max and mean over E1-E5 as
above):

| N | Vm L2 | Vm Linf | order L2 | order Linf | ECG max | order | ECG mean | order |
|---|---|---|---|---|---|---|---|---|
| 10 | 1.3280e-03 | 4.3786e-03 | | | 2.1139e-03 | | 1.1463e-03 | |
| 20 | 3.5025e-04 | 1.1585e-03 | 1.92 | 1.92 | 6.2174e-04 | 1.77 | 3.3928e-04 | 1.76 |
| 40 | 8.8865e-05 | 2.9327e-04 | 1.98 | 1.98 | 1.6761e-04 | 1.89 | 9.1642e-05 | 1.89 |
| 80 | 2.2284e-05 | 7.3449e-05 | 2.00 | 2.00 | 4.3431e-05 | 1.95 | 2.3761e-05 | 1.95 |

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
