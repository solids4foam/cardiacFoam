# cartesianConvergence - eikonalECG

Spatial convergence on Cartesian hexahedra (`Allrun`, `-mesh hex`) with the
case defaults (rotated conductivity, advection-diffusion formulation,
`leastSquares`): the activation time in 1D, 2D and 3D, and the 3D
eikonal-driven ECG at E1-E5 against its reference integral.

## Ladder

| N | hex cells (1D / 2D / 3D) |
|---|---|
| 10 | 10 / 100 / 1000 |
| 20 | 20 / 400 / 8000 |
| 40 | 40 / 1600 / 64000 |
| 80 | 80 / 6400 / 512000 |

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`hex_<d>D_N<n>`, the names that `../summarize.py` reads. From this folder:

```bash
mkdir -p results
for n in 10 20 40 80
do
    rsync -a --exclude setup ../../../ results/hex_1D_N$n/
    results/hex_1D_N$n/Allrun -dim 1 -N $n
    for d in 2 3
    do
        rsync -a --exclude setup ../../../ results/hex_${d}D_N$n/
        results/hex_${d}D_N$n/Allrun parallel -dim $d -N $n
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_hex_convergence.json`.

## Measured behaviour

Results of the last run of this study (`results/`).

Activation time `psi` (outer iterations to the residual tolerance `1e-8`):

| dim | N | L2 | Linf | order L2 | order Linf | outer iterations |
|---|---|---|---|---|---|---|
| 1D | 10 | 1.5542e-03 | 4.3264e-03 | | | 30 |
| 1D | 20 | 2.7029e-04 | 9.5512e-04 | 2.52 | 2.18 | 83 |
| 1D | 40 | 5.5578e-05 | 2.2463e-04 | 2.28 | 2.09 | 269 |
| 1D | 80 | 1.2859e-05 | 5.4459e-05 | 2.11 | 2.04 | 940 |
| 2D | 10 | 3.3714e-03 | 1.2521e-02 | | | 35 |
| 2D | 20 | 6.2456e-04 | 2.5585e-03 | 2.43 | 2.29 | 99 |
| 2D | 40 | 1.2946e-04 | 5.6566e-04 | 2.27 | 2.18 | 328 |
| 2D | 80 | 2.8569e-05 | 1.3197e-04 | 2.18 | 2.10 | 1162 |
| 3D | 10 | 3.9439e-03 | 2.0407e-02 | | | 34 |
| 3D | 20 | 7.3250e-04 | 4.1344e-03 | 2.43 | 2.30 | 95 |
| 3D | 40 | 1.5993e-04 | 9.5348e-04 | 2.20 | 2.12 | 318 |
| 3D | 80 | 3.9015e-05 | 2.3094e-04 | 2.04 | 2.05 | 1128 |

The activation time approaches second order from above in every dimension.
The reconstruction of the exact gradient (`gradReconActivationTime`) converges
at 1.92, 1.96, 1.98 in 3D, its boundary part at about 2.4.

3D ECG, E1-E5 (largest and mean over the electrodes of each electrode's Linf
error over time):

| N | max | order | mean | order |
|---|---|---|---|---|
| 10 | 4.5212e-01 | | 3.0139e-01 | |
| 20 | 4.9196e-03 | 6.52 | 1.9839e-03 | 7.25 |
| 40 | 3.6381e-05 | 7.08 | 2.6979e-05 | 6.20 |
| 80 | 7.8317e-06 | 2.22 | 6.2339e-06 | 2.11 |

The template `U(t - psi)` completes about nineteen periods across the cube
(`psi` runs from 1 to `e^3`), which `N = 10` and `N = 20` do not resolve: the
`N = 10` error is of the order of the signal, and the first two orders are
pre-asymptotic. The finest pair is near second order.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
