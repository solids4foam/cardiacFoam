# tetConvergence - bidomain/sphericalShell

Spatial convergence on gmsh tetrahedra (`Allrun -mesh tet`): sealed curved
walls approximated by planar faces, two unequal oblique full conductivity
tensors and skewed wall cells.

## Ladder

| N | lc | cells (three seeds) |
|---|---|---|
| 6 | 0.167 | 4430 - 4500 |
| 8 | 0.125 | 9815 - 9963 |
| 10 | 0.1 | 18059 - 18177 |
| 14 | 0.0714 | 46846 - 47129 |
| 20 | 0.05 | 137124 - 138569 |
| 28 | 0.0357 | 368100 - 368580 |
| 40 | 0.025 | 1071137 - 1075610 |

Meshes at different `N` are generated independently (gmsh 4.15.2). Every
level runs with `-seed 1`, `-seed 2` and `-seed 3` for both wall gradients;
the T = 0.2 tables give the mean over the seeds. With `-seed 1`: `T = 2` at
`N = 10, 14, 20, 28`, half the default time step at `N = 14, 28`, the `phiE`
tolerance `1e-6` at `N = 20, 28`, and the end-time controls at
`N = 14, 20, 28`. Every run without `-dt` uses the `Allrun` default time step,
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which scales as `h^2` and ends
exactly at `endTime`. The verifier writes `postProcessing/3D_tet_<N>_cells.dat`
(rows and header lines as in [hexConvergence](../hexConvergence/README.md)).

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<sealed|stock>_tet_N<n>_s<seed>[_T2|_dthalf|_tol1e-6]`, the names that
`../summarize.py` reads: the wall gradient (`-wallGradient`), `N`, the gmsh
seed, the long window (`-endTime 2`), half the default time step at
`T = 0.2` (`-dt`) or the `phiE` tolerance (`-phiTolerance`). The end-time
controls go to `results/endTimeControls/sealed_tet_N<n>_T<T>` (end time `T`)
and `sealed_tet_N<n>_steps<k>` (`k` default time steps). From this folder:

```bash
mkdir -p results/endTimeControls
for g in sealed stock
do
    for n in 6 8 10 14 20 28 40
    do
        for s in 1 2 3
        do
            rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_s$s/
            results/${g}_tet_N${n}_s$s/Allrun parallel -mesh tet -N $n -seed $s -wallGradient $g
        done
    done
    for n in 10 14 20 28
    do
        rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_s1_T2/
        results/${g}_tet_N${n}_s1_T2/Allrun parallel -mesh tet -N $n -seed 1 -wallGradient $g -endTime 2
    done
    for n in 14 28
    do
        dt=$(awk -v N=$n 'BEGIN { h = 0.00892857*(10/N)^2; k = int(0.2/h); if (k*h < 0.2) k++; printf "%.17g", 0.1/k }')
        rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_s1_dthalf/
        results/${g}_tet_N${n}_s1_dthalf/Allrun parallel -mesh tet -N $n -seed 1 -wallGradient $g -dt $dt
    done
    for n in 20 28
    do
        rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_s1_tol1e-6/
        results/${g}_tet_N${n}_s1_tol1e-6/Allrun parallel -mesh tet -N $n -seed 1 -wallGradient $g -phiTolerance 1e-6
    done
done
for n in 14 20 28
do
    for T in 0.05 0.1
    do
        r=results/endTimeControls/sealed_tet_N${n}_T$T
        rsync -a --exclude setup ../../../ $r/
        $r/Allrun parallel -mesh tet -N $n -seed 1 -endTime $T
    done
    dt=$(awk -v N=$n 'BEGIN { h = 0.00892857*(10/N)^2; k = int(0.2/h); if (k*h < 0.2) k++; printf "%.17g", 0.2/k }')
    for k in 1 4
    do
        r=results/endTimeControls/sealed_tet_N${n}_steps$k
        rsync -a --exclude setup ../../../ $r/
        $r/Allrun parallel -mesh tet -N $n -seed 1 -dt $dt -endTime $(awk -v d=$dt -v k=$k 'BEGIN { printf "%.17g", d*k }')
    done
done
python3 ../summarize.py results
```

`-wallGradient stock` replaces the four `sealedWallLeastSquares` entries of
`system/fvSchemes` by `leastSquares`. `../summarize.py` prints every table
below except the end-time controls, which are read from the
`3D_tet_<N>_cells.dat` summaries of `results/endTimeControls/`.

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`. Errors are L2 norms against the exact fields, except the
Linf columns; `h` is the verifier's `# h`. `phiE_gauge` compares `phiE` after
removing the volume-weighted mean error. The boundary rows are normalised by
the total volume, so their order is half an order above the RMS order over the
wall band. Orders are least-squares fits of `log(error)` against `log(h)` over
the three finest levels (every seed at T = 0.2), with their standard errors.

### `sealedWallLeastSquares` (as shipped), T = 0.2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 6 | 0.0934 | 7.064e-03 | 5.028e-02 | 2.247e-01 | 5.089e-01 | 4.866e-03 | 3.177e-01 |
| 8 | 0.0718 | 5.146e-03 | 3.954e-02 | 1.048e-01 | 2.480e-01 | 3.419e-03 | 2.515e-01 |
| 10 | 0.0587 | 4.020e-03 | 3.535e-02 | 6.814e-02 | 1.752e-01 | 2.541e-03 | 2.030e-01 |
| 14 | 0.0427 | 2.340e-03 | 2.249e-02 | 3.380e-02 | 8.952e-02 | 1.436e-03 | 1.446e-01 |
| 20 | 0.0298 | 1.253e-03 | 1.590e-02 | 1.586e-02 | 4.388e-02 | 7.234e-04 | 9.626e-02 |
| 28 | 0.0215 | 6.833e-04 | 1.170e-02 | 7.990e-03 | 2.246e-02 | 3.796e-04 | 6.577e-02 |
| 40 | 0.0151 | 3.432e-04 | 7.818e-03 | 3.887e-03 | 1.122e-02 | 1.859e-04 | 4.217e-02 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 20, 28, 40 | 1.89 ± 0.01 | 1.03 ± 0.14 | 2.05 ± 0.01 | 1.99 ± 0.02 | 1.99 ± 0.01 | 1.21 ± 0.01 |

### `sealedWallLeastSquares`, T = 2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 10 | 0.0586 | 6.129e-02 | 3.995e-01 | 1.124e-01 | 2.934e-01 | 3.863e-02 | 4.770e-01 |
| 14 | 0.0428 | 2.902e-02 | 1.751e-01 | 5.629e-02 | 1.543e-01 | 1.579e-02 | 2.783e-01 |
| 20 | 0.0299 | 1.390e-02 | 7.571e-02 | 2.628e-02 | 7.735e-02 | 6.316e-03 | 1.698e-01 |
| 28 | 0.0215 | 7.094e-03 | 4.326e-02 | 1.322e-02 | 3.820e-02 | 2.703e-03 | 1.098e-01 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 14, 20, 28 | 2.05 ± 0.00 | 2.04 ± 0.19 | 2.11 ± 0.01 | 2.03 ± 0.06 | 2.57 ± 0.01 | 1.35 ± 0.02 |

### Stock `leastSquares` (control), T = 0.2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 6 | 0.0934 | 6.933e-03 | 3.869e-02 | 2.125e-01 | 4.753e-01 | 4.486e-03 | 2.453e-01 |
| 8 | 0.0718 | 5.056e-03 | 3.410e-02 | 9.535e-02 | 2.267e-01 | 3.198e-03 | 2.090e-01 |
| 10 | 0.0587 | 4.064e-03 | 3.319e-02 | 6.139e-02 | 1.573e-01 | 2.561e-03 | 1.916e-01 |
| 14 | 0.0427 | 2.580e-03 | 2.491e-02 | 2.995e-02 | 7.754e-02 | 1.722e-03 | 1.655e-01 |
| 20 | 0.0298 | 1.562e-03 | 2.044e-02 | 1.470e-02 | 3.683e-02 | 1.085e-03 | 1.425e-01 |
| 28 | 0.0215 | 9.978e-04 | 1.762e-02 | 8.135e-03 | 2.581e-02 | 7.242e-04 | 1.239e-01 |
| 40 | 0.0151 | 6.255e-04 | 1.334e-02 | 4.705e-03 | 1.771e-02 | 4.666e-04 | 1.063e-01 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 20, 28, 40 | 1.34 ± 0.01 | 0.61 ± 0.14 | 1.66 ± 0.03 | 1.07 ± 0.02 | 1.23 ± 0.01 | 0.43 ± 0.00 |

### Stock `leastSquares` (control), T = 2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 10 | 0.0586 | 6.043e-02 | 2.893e-01 | 1.013e-01 | 2.626e-01 | 3.639e-02 | 3.681e-01 |
| 14 | 0.0428 | 3.010e-02 | 1.840e-01 | 5.011e-02 | 1.366e-01 | 1.657e-02 | 2.804e-01 |
| 20 | 0.0299 | 1.605e-02 | 9.840e-02 | 2.446e-02 | 6.405e-02 | 7.796e-03 | 2.316e-01 |
| 28 | 0.0215 | 9.804e-03 | 6.542e-02 | 1.348e-02 | 4.556e-02 | 4.160e-03 | 1.985e-01 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 14, 20, 28 | 1.63 ± 0.07 | 1.51 ± 0.15 | 1.91 ± 0.06 | 1.61 ± 0.31 | 2.01 ± 0.06 | 0.50 ± 0.02 |

With `sealedWallLeastSquares` `phiE` is second order in L2 and Linf and `Vm`
reaches 1.89 in L2 at T = 0.2 and 2.05 at T = 2; with stock `leastSquares`
the `Vm` L2 order is 1.34 and the wall-cell gradient order 0.43.

### Seeds

Spread of the L2 errors over the three seeds at one `N`
(`max/min - 1`):

| N | sealed: Vm L2 | sealed: phiE_gauge L2 | stock: Vm L2 | stock: phiE_gauge L2 |
|---|---|---|---|---|
| 6 | 0.60 % | 8.01 % | 2.86 % | 8.68 % |
| 8 | 3.28 % | 2.55 % | 1.68 % | 2.27 % |
| 10 | 1.57 % | 5.39 % | 2.02 % | 6.06 % |
| 14 | 0.87 % | 3.37 % | 1.82 % | 2.63 % |
| 20 | 1.07 % | 2.35 % | 0.13 % | 1.86 % |
| 28 | 0.81 % | 1.92 % | 1.82 % | 1.58 % |
| 40 | 0.90 % | 0.58 % | 1.40 % | 0.46 % |

### Time step and `phiE` tolerance

Change of the T = 0.2 L2 errors under half the default time step and under
`-phiTolerance 1e-6` (`-seed 1`):

| gradient | N | dt/2: Vm L2 | dt/2: phiE_gauge L2 | tolerance 1e-6: Vm L2 | tolerance 1e-6: phiE_gauge L2 |
|---|---|---|---|---|---|
| `sealedWallLeastSquares` | 14 | +0.23 % | +0.00 % |  |  |
| `sealedWallLeastSquares` | 20 |  |  | -0.000 % | -0.006 % |
| `sealedWallLeastSquares` | 28 | +0.02 % | -0.03 % | +0.000 % | -0.036 % |
| stock `leastSquares` | 14 | +0.23 % | -0.05 % |  |  |
| stock `leastSquares` | 20 |  |  | -0.000 % | -0.007 % |
| stock `leastSquares` | 28 | +0.02 % | -0.11 % | -0.000 % | +0.001 % |

Both change the errors by at most 0.23 %, below the seed spread at the same
`N`.

### End time

`sealedWallLeastSquares`, `-seed 1`, `Vm` L2 and orders fitted over
`N = 14, 20, 28` at the end of each window; the one- and four-step windows
end after one and four default time steps, which shrink as `h^2`, so their
levels stop at different times. The `0.2` and `2` rows are the seed-1 runs
above.

| window | end time | Vm L2, N = 14 | Vm L2, N = 20 | Vm L2, N = 28 | order Vm L2 | order Vm Linf | order phiE L2 | order gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|---|
| 1 step | 0.004545, 0.002222, 0.001136 | 2.256e-04 | 1.075e-04 | 5.284e-05 | 2.11 | 1.11 | 2.18 | 1.23 |
| 4 steps | 0.01818, 0.008889, 0.004545 | 6.949e-04 | 3.294e-04 | 1.630e-04 | 2.11 | 1.23 | 2.12 | 1.23 |
| 0.05 | 0.05, 0.05, 0.05 | 1.229e-03 | 7.713e-04 | 4.652e-04 | 1.41 | 0.69 | 2.10 | 1.15 |
| 0.1 | 0.1, 0.1, 0.1 | 1.697e-03 | 9.881e-04 | 5.596e-04 | 1.61 | 0.84 | 2.10 | 1.14 |
| 0.2 | 0.2, 0.2, 0.2 | 2.327e-03 | 1.260e-03 | 6.799e-04 | 1.79 | 1.09 | 2.10 | 1.14 |
| 2 | 2, 2, 2 | 2.902e-02 | 1.390e-02 | 7.094e-03 | 2.05 | 2.04 | 2.11 | 1.35 |

The `phiE` L2 order is 2.10-2.18 and the wall-cell gradient order 1.14-1.35
in every window. The `Vm` L2 order is 2.11 after one and four steps, 1.41 at
`T = 0.05`, 1.79 at `T = 0.2` and 2.05 at `T = 2`.

### Linear solves

Every linear solve reaches its tolerance below the explicit `maxIter 5000`:

| runs | field | max iterations | max final residual |
|---|---|---|---|
| sealed tet T = 0.2 | Vm | 14 | 1.00e-15 |
| sealed tet T = 0.2 | phiE | 836 | 1.00e-15 |
| sealed tet T2 | Vm | 14 | 9.99e-16 |
| sealed tet T2 | phiE | 601 | 1.00e-15 |
| sealed tet dthalf | Vm | 10 | 9.92e-16 |
| sealed tet dthalf | phiE | 595 | 1.00e-15 |
| sealed tet tol1e-6 | Vm | 13 | 9.94e-16 |
| sealed tet tol1e-6 | phiE | 331 | 9.99e-07 |
| stock tet T = 0.2 | Vm | 14 | 1.00e-15 |
| stock tet T = 0.2 | phiE | 849 | 1.00e-15 |
| stock tet T2 | Vm | 13 | 1.00e-15 |
| stock tet T2 | phiE | 605 | 1.00e-15 |
| stock tet dthalf | Vm | 10 | 9.99e-16 |
| stock tet dthalf | phiE | 599 | 1.00e-15 |
| stock tet tol1e-6 | Vm | 14 | 9.43e-16 |
| stock tet tol1e-6 | phiE | 317 | 1.00e-06 |

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
