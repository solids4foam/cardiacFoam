# hexConvergence - monodomainPseudoECG/sphericalShell

Spatial convergence on the cubed sphere (`Allrun -mesh hex`, default): sealed
curved walls approximated by planar faces, an oblique full conductivity tensor
and nearly orthogonal wall cells.

## Ladder

| N | k | cells |
|---|---|---|
| 8 | 2 | 864 |
| 12 | 3 | 2916 |
| 16 | 4 | 6912 |
| 24 | 6 | 23328 |
| 32 | 8 | 55296 |
| 40 | 10 | 108000 |
| 48 | 12 | 186624 |
| 64 | 16 | 442368 |

`N` is a multiple of 4 and `k = N/4`. The verifier writes
`postProcessing/3D_hex_<N>_cells.dat`: rows `Vm`, `u1`, `u2`, the bulk,
boundary and gradient rows, and the header lines `# h`, `# time` and
`# cellsPerDirection`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<sealed|stock>_hex_N<n>_s<seed>[_T2|_dthalf]`, the names that
`../summarize.py` reads. From this folder, for each `N`:

```bash
mkdir -p results
rsync -a --exclude setup ../../../ results/sealed_hex_N<n>_s1/
results/sealed_hex_N<n>_s1/Allrun parallel -mesh hex -N <n>
```

The default end time is `T = 0.2` and the time step is
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which scales as `h^2`.

- Stock control: copies named `stock_...`, with `-wallGradient stock`, which
  replaces the two `sealedWallLeastSquares` entries of `system/fvSchemes`
  by `leastSquares`.
- Long window: copies ending in `_T2`, with `-endTime 2`.
- Time-step control: copies ending in `_dthalf`, with `-dt` set to half the
  default time step at `T = 0.2`.

Every table below is printed by

```bash
python3 ../summarize.py results
```

## Measured behaviour

Errors are L2 norms against the exact fields, except the `Vm` Linf column;
`h` is the verifier's `# h`. The boundary rows are normalised by the total
volume, so their order is half an order above the RMS order over the wall
band. Orders are least-squares fits of `log(error)` against `log(h)` over the
three finest levels, with their standard errors.

### `sealedWallLeastSquares` (as shipped), T = 0.2

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 8 | 0.1603 | 5.148e-03 | 2.659e-02 | 4.682e-03 | 3.401e-01 | 3.216e-01 |
| 12 | 0.1074 | 2.912e-03 | 1.737e-02 | 2.558e-03 | 2.550e-01 | 2.438e-01 |
| 16 | 0.0807 | 1.775e-03 | 1.238e-02 | 1.535e-03 | 1.904e-01 | 1.826e-01 |
| 24 | 0.0539 | 8.113e-04 | 7.005e-03 | 6.716e-04 | 1.180e-01 | 1.132e-01 |
| 32 | 0.0404 | 4.545e-04 | 4.464e-03 | 3.566e-04 | 8.142e-02 | 7.813e-02 |
| 40 | 0.0324 | 2.896e-04 | 3.079e-03 | 2.154e-04 | 6.025e-02 | 5.791e-02 |
| 48 | 0.0270 | 2.006e-04 | 2.357e-03 | 1.419e-04 | 4.680e-02 | 4.508e-02 |
| 64 | 0.0202 | 1.127e-04 | 1.595e-03 | 7.297e-05 | 3.111e-02 | 3.010e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 40, 48, 64 | 2.01 ± 0.00 | 1.40 ± 0.03 | 2.31 ± 0.01 | 1.41 ± 0.01 | 1.39 ± 0.01 |

`Vm` is second order in L2. The wall-cell gradient is first order (the
boundary rows less half an order); `gradReconVm_boundary`, the scheme applied
to the exact field, shows that this order belongs to the reconstruction, which
takes no data from the sealed faces.

### `sealedWallLeastSquares`, T = 2

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 16 | 0.0807 | 4.702e-02 | 1.984e-01 | 3.134e-02 | 4.062e-01 | 2.887e-01 |
| 24 | 0.0539 | 2.191e-02 | 8.845e-02 | 1.203e-02 | 2.237e-01 | 1.789e-01 |
| 32 | 0.0404 | 1.265e-02 | 5.196e-02 | 6.045e-03 | 1.467e-01 | 1.235e-01 |
| 48 | 0.0270 | 5.773e-03 | 2.382e-02 | 2.267e-03 | 8.027e-02 | 7.127e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 24, 32, 48 | 1.93 ± 0.01 | 1.90 ± 0.02 | 2.41 ± 0.01 | 1.48 ± 0.01 | 1.33 ± 0.02 |

### Stock `leastSquares` (control), T = 0.2

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 8 | 0.1603 | 4.468e-03 | 2.699e-02 | 3.940e-03 | 2.402e-01 | 1.382e-01 |
| 12 | 0.1074 | 2.475e-03 | 1.750e-02 | 2.052e-03 | 1.846e-01 | 7.605e-02 |
| 16 | 0.0807 | 1.654e-03 | 1.132e-02 | 1.387e-03 | 1.624e-01 | 4.491e-02 |
| 24 | 0.0539 | 1.008e-03 | 6.762e-03 | 8.793e-04 | 1.385e-01 | 1.928e-02 |
| 32 | 0.0404 | 7.452e-04 | 6.242e-03 | 6.526e-04 | 1.235e-01 | 1.009e-02 |
| 40 | 0.0324 | 5.981e-04 | 5.725e-03 | 5.153e-04 | 1.125e-01 | 6.002e-03 |
| 48 | 0.0270 | 5.016e-04 | 5.345e-03 | 4.214e-04 | 1.039e-01 | 3.896e-03 |
| 64 | 0.0202 | 3.807e-04 | 4.586e-03 | 3.016e-04 | 9.128e-02 | 1.950e-03 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 40, 48, 64 | 0.96 ± 0.00 | 0.48 ± 0.04 | 1.14 ± 0.02 | 0.45 ± 0.00 | 2.39 ± 0.01 |

`Vm` is first order: the zero-gradient wall value enters the wall-cell
gradient, whose error does not vanish with `h` (the boundary rows less half an
order). With exact wall values (`gradReconVm_boundary`) the same scheme
converges, so the error comes from the wall value, not from the stencil.

### Stock `leastSquares` (control), T = 2

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 16 | 0.0807 | 4.533e-02 | 1.920e-01 | 3.088e-02 | 3.201e-01 | 7.101e-02 |
| 24 | 0.0539 | 2.552e-02 | 1.024e-01 | 1.492e-02 | 2.522e-01 | 3.049e-02 |
| 32 | 0.0404 | 1.816e-02 | 8.363e-02 | 9.624e-03 | 2.168e-01 | 1.596e-02 |
| 48 | 0.0270 | 1.191e-02 | 6.320e-02 | 5.485e-03 | 1.757e-01 | 6.159e-03 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 24, 32, 48 | 1.10 ± 0.04 | 0.70 ± 0.00 | 1.44 ± 0.04 | 0.52 ± 0.00 | 2.31 ± 0.03 |

### Time step

Halving the time step at T = 0.2 changes the `Vm` L2 error by -0.14 %
(`N = 32`) and -0.24 % (`N = 64`) with `sealedWallLeastSquares`, and by
-0.21 % and -0.17 % with the stock gradient: the errors above are spatial.

### `-trace conormal` (stock `leastSquares`)

Runs with `-trace conormal -wallGradient stock`; `summarize.py` does not read
them, and the orders are between successive levels.

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 8 | 0.1603 | 4.139e-03 | 2.601e-02 | 3.574e-03 | 2.774e-01 | 1.382e-01 |
| 12 | 0.1074 | 2.007e-03 | 1.587e-02 | 1.454e-03 | 1.386e-01 | 7.605e-02 |
| 16 | 0.0807 | 1.156e-03 | 9.221e-03 | 7.241e-04 | 8.106e-02 | 4.491e-02 |
| 24 | 0.0539 | 5.224e-04 | 4.200e-03 | 2.671e-04 | 3.818e-02 | 1.928e-02 |
| 32 | 0.0404 | 2.972e-04 | 2.639e-03 | 1.356e-04 | 2.314e-02 | 1.009e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 8 -> 12 | 1.81 | 1.24 | 2.25 | 1.73 | 1.49 |
| 12 -> 16 | 1.93 | 1.90 | 2.44 | 1.88 | 1.84 |
| 16 -> 24 | 1.97 | 1.95 | 2.47 | 1.86 | 2.09 |
| 24 -> 32 | 1.96 | 1.62 | 2.36 | 1.74 | 2.25 |

## Serial and decomposed runs

With `sealedWallLeastSquares`, a serial run and a four-rank scotch run
(`parallel`) at `N = 16` agree to a maximum relative difference of 3.1e-12
over all rows.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
