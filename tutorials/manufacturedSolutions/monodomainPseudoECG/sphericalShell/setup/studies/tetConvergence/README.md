# tetConvergence - monodomainPseudoECG/sphericalShell

Spatial convergence on gmsh tetrahedra (`Allrun -mesh tet`): sealed curved
walls approximated by planar faces, an oblique full conductivity tensor and
skewed wall cells.

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

Meshes at different `N` are generated independently. Every level runs with
`-seed 1`, `-seed 2` and `-seed 3`; the tables give the mean over the seeds.
The verifier writes `postProcessing/3D_tet_<N>_cells.dat` (rows and header
lines as in [hexConvergence](../hexConvergence/README.md)).

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<sealed|stock>_tet_N<n>_s<seed>[_T2|_dthalf]`, the names that
`../summarize.py` reads. From this folder, for each `N` and seed `s`:

```bash
mkdir -p results
rsync -a --exclude setup ../../../ results/sealed_tet_N<n>_s<s>/
results/sealed_tet_N<n>_s<s>/Allrun parallel -mesh tet -N <n> -seed <s>
```

The default end time is `T = 0.2` and the time step is
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which scales as `h^2`.

- Stock control: copies named `stock_...`, with `-wallGradient stock`, which
  replaces the two `sealedWallLeastSquares` entries of `system/fvSchemes`
  by `leastSquares`.
- Long window: copies ending in `_T2`, with `-endTime 2` (`-seed 1`).
- Time-step control: copies ending in `_dthalf`, with `-dt` set to half the
  default time step at `T = 0.2`.

Every table below is printed by

```bash
python3 ../summarize.py results
```

## Measured behaviour

Errors are L2 norms against the exact fields, except the `Vm` Linf column;
`h` is the verifier's `# h`. At T = 0.2 the values are means over the three
seeds; the orders are least-squares fits of `log(error)` against `log(h)` over
every seed of the three finest levels, with their standard errors. The
boundary rows are normalised by the total volume, so their order is half an
order above the RMS order over the wall band.

### `sealedWallLeastSquares` (as shipped), T = 0.2

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 6 | 0.0934 | 6.005e-03 | 4.366e-02 | 4.634e-03 | 2.618e-01 | 2.474e-01 |
| 8 | 0.0718 | 4.430e-03 | 4.233e-02 | 3.344e-03 | 2.195e-01 | 2.030e-01 |
| 10 | 0.0587 | 3.513e-03 | 3.695e-02 | 2.532e-03 | 1.851e-01 | 1.694e-01 |
| 14 | 0.0427 | 2.100e-03 | 2.466e-02 | 1.465e-03 | 1.386e-01 | 1.226e-01 |
| 20 | 0.0298 | 1.144e-03 | 2.060e-02 | 7.575e-04 | 9.632e-02 | 8.215e-02 |
| 28 | 0.0215 | 6.274e-04 | 1.377e-02 | 3.998e-04 | 6.626e-02 | 5.455e-02 |
| 40 | 0.0151 | 3.174e-04 | 8.536e-03 | 1.975e-04 | 4.240e-02 | 3.392e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 20, 28, 40 | 1.87 ± 0.01 | 1.27 ± 0.14 | 1.96 ± 0.01 | 1.20 ± 0.01 | 1.29 ± 0.01 |

Over the seeds at one `N` the `Vm` L2 error varies by at most 1.4 %
(1.139e-03 to 1.154e-03 at `N = 20`; 3.165e-04 to 3.190e-04 at `N = 40`).

### `sealedWallLeastSquares`, T = 2

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 10 | 0.0586 | 5.009e-02 | 3.285e-01 | 3.350e-02 | 4.182e-01 | 2.677e-01 |
| 14 | 0.0428 | 2.368e-02 | 1.462e-01 | 1.321e-02 | 2.588e-01 | 1.939e-01 |
| 20 | 0.0299 | 1.167e-02 | 7.590e-02 | 5.299e-03 | 1.678e-01 | 1.306e-01 |
| 28 | 0.0215 | 6.140e-03 | 4.371e-02 | 2.332e-03 | 1.099e-01 | 8.608e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 14, 20, 28 | 1.96 ± 0.01 | 1.76 ± 0.04 | 2.52 ± 0.02 | 1.25 ± 0.02 | 1.18 ± 0.05 |

### Stock `leastSquares` (control), T = 0.2

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 6 | 0.0934 | 6.111e-03 | 4.217e-02 | 4.530e-03 | 2.057e-01 | 1.543e-01 |
| 8 | 0.0718 | 4.364e-03 | 3.679e-02 | 3.153e-03 | 1.743e-01 | 1.143e-01 |
| 10 | 0.0587 | 3.528e-03 | 3.387e-02 | 2.483e-03 | 1.606e-01 | 9.547e-02 |
| 14 | 0.0427 | 2.306e-03 | 2.475e-02 | 1.678e-03 | 1.410e-01 | 7.222e-02 |
| 20 | 0.0298 | 1.462e-03 | 2.004e-02 | 1.096e-03 | 1.233e-01 | 5.534e-02 |
| 28 | 0.0215 | 9.721e-04 | 1.688e-02 | 7.544e-04 | 1.082e-01 | 4.501e-02 |
| 40 | 0.0151 | 6.355e-04 | 1.354e-02 | 5.020e-04 | 9.326e-02 | 3.662e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 20, 28, 40 | 1.22 ± 0.01 | 0.57 ± 0.05 | 1.14 ± 0.01 | 0.41 ± 0.00 | 0.60 ± 0.01 |

### Stock `leastSquares` (control), T = 2

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 10 | 0.0586 | 4.791e-02 | 3.094e-01 | 3.054e-02 | 3.051e-01 | 1.516e-01 |
| 14 | 0.0428 | 2.490e-02 | 1.667e-01 | 1.401e-02 | 2.391e-01 | 1.148e-01 |
| 20 | 0.0299 | 1.517e-02 | 1.027e-01 | 7.392e-03 | 2.020e-01 | 8.721e-02 |
| 28 | 0.0215 | 1.001e-02 | 7.473e-02 | 4.226e-03 | 1.740e-01 | 7.109e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 14, 20, 28 | 1.33 ± 0.04 | 1.17 ± 0.11 | 1.75 ± 0.03 | 0.46 ± 0.01 | 0.70 ± 0.04 |

### Time step

Halving the time step at T = 0.2 (`-seed 1`) changes the `Vm` L2 error by
+0.11 % (`N = 14`) and -0.14 % (`N = 40`) with `sealedWallLeastSquares`, and
by +0.10 % and -0.13 % with the stock gradient.

### `-trace conormal` (stock `leastSquares`)

Runs with `-trace conormal -wallGradient stock`; `summarize.py` does not read
them, and the orders are between successive levels.

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 6 | 0.0934 | 1.533e-02 | 1.435e+00 | 1.362e-02 | 8.401e+00 | 1.543e-01 |
| 8 | 0.0718 | 5.793e-03 | 9.625e-02 | 4.470e-03 | 2.252e-01 | 1.143e-01 |
| 10 | 0.0587 | 4.496e-03 | 7.338e-02 | 3.365e-03 | 1.816e-01 | 9.547e-02 |
| 14 | 0.0427 | 5.417e-03 | 8.927e-01 | 2.492e-03 | 8.663e+00 | 7.222e-02 |

These means include three meshes with a local wall-gradient blow-up, so the
orders between levels are not meaningful (they range from -12.17 to +13.77
for `gradVm_boundary`). Per seed:

| run | Vm L2 | Vm Linf | gradVm_boundary L2 | gradVm_boundary Linf |
|---|---|---|---|---|
| N 6, seed 1 | 8.36e-03 | 9.76e-02 | 3.42e-01 | 1.29e+01 |
| N 6, seed 2 | 8.05e-03 | 7.19e-02 | 3.03e-01 | 6.85e+00 |
| N 6, seed 3 | 2.96e-02 | 4.14e+00 | 2.46e+01 | 3.75e+03 |
| N 8, seed 1 | 6.04e-03 | 7.48e-02 | 2.25e-01 | 4.10e+00 |
| N 8, seed 2 | 5.65e-03 | 9.15e-02 | 2.17e-01 | 3.45e+00 |
| N 8, seed 3 | 5.68e-03 | 1.22e-01 | 2.34e-01 | 1.12e+01 |
| N 10, seed 1 | 4.48e-03 | 7.30e-02 | 1.71e-01 | 3.88e+00 |
| N 10, seed 2 | 4.49e-03 | 9.04e-02 | 2.02e-01 | 1.64e+01 |
| N 10, seed 3 | 4.52e-03 | 5.67e-02 | 1.71e-01 | 3.44e+00 |
| N 14, seed 1 | 2.54e-03 | 4.51e-02 | 1.09e-01 | 3.22e+00 |
| N 14, seed 2 | 4.66e-03 | 7.62e-01 | 7.89e+00 | 2.26e+03 |
| N 14, seed 3 | 9.05e-03 | 1.87e+00 | 1.80e+01 | 5.38e+03 |

Three of the twelve meshes (`N = 6` seed 3, `N = 14` seeds 2 and 3) show a
local wall-gradient blow-up, with `gradVm_boundary` Linf of order 1e3; on the
other nine `gradVm_boundary` Linf stays between 3.2 and 16.4.
`gradReconVm_boundary`, computed from the exact field, is unaffected.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
