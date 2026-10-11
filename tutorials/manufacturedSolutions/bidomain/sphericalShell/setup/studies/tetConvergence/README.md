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

Meshes at different `N` are generated independently. Every level runs with
`-seed 1`, `-seed 2` and `-seed 3`; the T = 0.2 tables give the mean over the
seeds. The verifier writes `postProcessing/3D_tet_<N>_cells.dat` (rows and
header lines as in [hexConvergence](../hexConvergence/README.md)).

## Execution

In a copy of the case, for each `N` and `s`:

```bash
./Allrun parallel -mesh tet -N <n> -seed <s>
```

The default end time is `T = 0.2` and the time step is
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`.

The case selects `sealedWallLeastSquares` for `grad(Vm)`, `grad(Vm.component(0))`,
`grad(phiE)` and `grad(phiE.component(0))` in `system/fvSchemes`. The stock wall
gradient is the control: in the run copy, delete those four lines, so every
gradient is `leastSquares`.

### Long window

The long window runs to `T = 2` (`-seed 1`); the time-step control repeats
the short window with the default time step halved:

```bash
./Allrun parallel -mesh tet -N <n> -seed 1 -endTime 2
./Allrun parallel -mesh tet -N <n> -seed 1 -dt <half the default>
```

## Measured behaviour

Errors are L2 norms against the exact fields, except the Linf columns; `h` is
the verifier's `# h`. `phiE_gauge` compares `phiE` after removing the
volume-weighted mean error. The boundary rows are normalised by the total
volume, so their order is half an order above the RMS order over the wall
band. Orders are least-squares fits of `log(error)` against `log(h)` over the
three finest levels (every seed at T = 0.2), with their standard errors.

### `sealedWallLeastSquares` (as shipped), T = 0.2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 6 | 0.0934 | 7.064e-03 | 5.028e-02 | 2.247e-01 | 5.089e-01 | 4.866e-03 | 3.177e-01 |
| 8 | 0.0718 | 5.146e-03 | 3.954e-02 | 1.048e-01 | 2.480e-01 | 3.419e-03 | 2.515e-01 |
| 10 | 0.0587 | 4.020e-03 | 3.535e-02 | 6.814e-02 | 1.752e-01 | 2.541e-03 | 2.030e-01 |
| 14 | 0.0427 | 2.340e-03 | 2.249e-02 | 3.380e-02 | 8.952e-02 | 1.436e-03 | 1.446e-01 |
| 20 | 0.0298 | 1.253e-03 | 1.590e-02 | 1.586e-02 | 4.388e-02 | 7.234e-04 | 9.626e-02 |
| 28 | 0.0215 | 6.833e-04 | 1.170e-02 | 7.990e-03 | 2.246e-02 | 3.796e-04 | 6.577e-02 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 14, 20, 28 | 1.79 ± 0.01 | 0.96 ± 0.16 | 2.10 ± 0.01 | 2.02 ± 0.02 | 1.94 ± 0.02 | 1.15 ± 0.01 |

Over the seeds at one `N` the `Vm` L2 error varies by at most 1.1 % and the
`phiE_gauge` L2 error by at most 2.4 % (`N = 20` and `28`).

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

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 14, 20, 28 | 1.39 ± 0.01 | 0.50 ± 0.15 | 1.90 ± 0.02 | 1.61 ± 0.11 | 1.26 ± 0.02 | 0.42 ± 0.00 |

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

### Time step

Halving the time step at `N = 14`, T = 0.2 (`-seed 1`), changes the `Vm` L2
error by +0.23 % and the `phiE_gauge` L2 error by less than 0.1 % with either
gradient.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
