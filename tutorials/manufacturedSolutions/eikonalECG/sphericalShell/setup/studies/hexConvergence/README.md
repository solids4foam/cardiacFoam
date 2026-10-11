# hexConvergence - eikonalECG/sphericalShell

Spatial convergence and outer-iteration behaviour on the cubed sphere
(`Allrun -mesh hex`, default): sealed curved walls approximated by planar
faces, an oblique full conductivity tensor and nearly orthogonal wall cells.

## Ladder

| N | k | cells |
|---|---|---|
| 8 | 2 | 864 |
| 12 | 3 | 2916 |
| 16 | 4 | 6912 |
| 24 | 6 | 23328 |
| 32 | 8 | 55296 |
| 48 | 12 | 186624 |

`N` is a multiple of 4 and `k = N/4`. The verifier writes
`postProcessing/3D_hex_<N>_cells.dat`: rows `activationTime`, the bulk,
boundary and gradient rows, and the header lines `# h`, `# time`,
`# cellsPerDirection` and `# pinnedCells`. The pinned ball (radius 0.15
around the stimulus point) holds 4, 11, 24, 80, 193 and 660 cells at
`N = 8, 12, 16, 24, 32, 48`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<sealed|stock>_hex_N<n>_s1[_ex][_exact|_best|_cap]`, the names that
`../summarize.py` reads: wall gradient, `N`, seed, explicit formulation, one
iteration from the exact field, stopped at the iteration of the smallest
outer residual, or stopped at an iteration cap. From this folder, for each `N`:

```bash
mkdir -p results
rsync -a --exclude setup ../../../ results/sealed_hex_N<n>_s1/
results/sealed_hex_N<n>_s1/Allrun -mesh hex -N <n>
```

- Stock gradient: copies named `stock_...`, with `-wallGradient stock`.
- Exact start: copies ending in `_exact`, with `-exactStart -iterations 1`.
- Explicit formulation: copies with `_ex_cap`, with `-formulation explicit
  -iterations 1500`.

Every table below is printed by

```bash
python3 ../summarize.py results
```

## Measured behaviour

Errors are L2 norms against the exact fields, except the `Linf` column; `h`
is the verifier's `# h`. `outer its` is the number of outer iterations,
`min residual` the smallest initial residual of the `activationTime` solves
with the iteration where it occurs, `last residual` that of the last solve.
Orders are least-squares fits of `log(error)` against `log(h)` over the three
finest levels, with their standard errors. `activationTime_boundary` is
normalised by the total volume, so its order is half an order above the RMS
order over the wall band.

With the case defaults (advection-diffusion formulation) both wall gradients
converge; the number of outer iterations grows roughly as `N^2`, and the
`N = 48` run stops at the 8000-iteration cap with the outer residual at
1.3e-8 (falling by about 2.8 per 1000 iterations). One iteration from the
exact field converges at a higher order than the converged field. With the
explicit formulation and `sealedWallLeastSquares` the outer iteration does
not converge: the outer residual passes a minimum and the activation time
drifts away from the pinned cells.

### sealed, hex, advectionDiffusion, converged

| N | h | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 | outer its | min residual | last residual |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 0.1603 | 2.818e+00 | 6.091e+00 | 2.067e+00 | 6.011e+00 | 3.837e+00 | 154 | 8.4e-09 @ 154 | 8.4e-09 |
| 12 | 0.1074 | 2.042e+00 | 4.921e+00 | 1.229e+00 | 4.581e+00 | 2.988e+00 | 303 | 9.1e-09 @ 303 | 9.1e-09 |
| 16 | 0.0807 | 1.604e+00 | 3.845e+00 | 8.353e-01 | 3.605e+00 | 2.428e+00 | 548 | 9.6e-09 @ 548 | 9.6e-09 |
| 24 | 0.0539 | 1.028e+00 | 2.406e+00 | 4.379e-01 | 2.332e+00 | 1.662e+00 | 1449 | 9.9e-09 @ 1449 | 9.9e-09 |
| 32 | 0.0404 | 7.022e-01 | 1.602e+00 | 2.597e-01 | 1.604e+00 | 1.207e+00 | 3069 | 9.9e-09 @ 3069 | 9.9e-09 |
| 48 | 0.0270 | 3.665e-01 | 8.024e-01 | 1.112e-01 | 8.949e-01 | 7.326e-01 | 8000 | 1.3e-08 @ 8000 | 1.3e-08 |

| N | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|
| 24, 32, 48 | 1.50 ± 0.08 | 1.59 ± 0.08 | 1.99 ± 0.08 | 1.39 ± 0.04 | 1.19 ± 0.03 |

### sealed, hex, advectionDiffusion, exact

| N | h | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 | outer its | min residual | last residual |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 0.1603 | 2.161e-01 | 1.051e+00 | 1.982e-01 | 4.157e+00 | 3.837e+00 | 1 | 2.6e-02 @ 1 | 2.6e-02 |
| 12 | 0.1074 | 8.804e-02 | 5.053e-01 | 7.738e-02 | 3.321e+00 | 2.988e+00 | 1 | 1.2e-02 @ 1 | 1.2e-02 |
| 16 | 0.0807 | 4.497e-02 | 3.077e-01 | 3.849e-02 | 2.706e+00 | 2.428e+00 | 1 | 6.7e-03 @ 1 | 6.7e-03 |
| 24 | 0.0539 | 1.641e-02 | 2.506e-01 | 1.396e-02 | 1.826e+00 | 1.662e+00 | 1 | 2.3e-03 @ 1 | 2.3e-03 |
| 32 | 0.0404 | 7.683e-03 | 1.784e-01 | 6.573e-03 | 1.305e+00 | 1.207e+00 | 1 | 9.4e-04 @ 1 | 9.4e-04 |
| 48 | 0.0270 | 2.487e-03 | 6.903e-02 | 2.137e-03 | 7.743e-01 | 7.326e-01 | 1 | 2.5e-04 @ 1 | 2.5e-04 |

| N | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|
| 24, 32, 48 | 2.73 ± 0.04 | 1.89 ± 0.32 | 2.72 ± 0.04 | 1.24 ± 0.03 | 1.19 ± 0.03 |

### sealed, hex, explicit, cap

| N | h | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 | outer its | min residual | last residual |
|---|---|---|---|---|---|---|---|---|---|
| 12 | 0.1074 | 4.708e+26 | 7.902e+26 | 2.804e+26 | 2.213e+26 | 2.988e+00 | 1500 | 1.5e-02 @ 19 | 7.1e-02 |
| 16 | 0.0807 | 2.672e+10 | 4.197e+10 | 1.379e+10 | 1.018e+10 | 2.428e+00 | 1500 | 4.8e-03 @ 45 | 4.0e-02 |

### stock, hex, advectionDiffusion, converged

| N | h | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 | outer its | min residual | last residual |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 0.1603 | 2.462e+00 | 6.280e+00 | 1.660e+00 | 8.664e+00 | 2.004e+00 | 173 | 8.7e-09 @ 173 | 8.7e-09 |
| 12 | 0.1074 | 1.712e+00 | 4.491e+00 | 9.432e-01 | 5.507e+00 | 1.050e+00 | 338 | 9.2e-09 @ 338 | 9.2e-09 |
| 16 | 0.0807 | 1.277e+00 | 3.388e+00 | 6.105e-01 | 4.064e+00 | 6.194e-01 | 607 | 9.6e-09 @ 607 | 9.6e-09 |
| 24 | 0.0539 | 8.128e-01 | 2.181e+00 | 3.184e-01 | 2.794e+00 | 2.752e-01 | 1631 | 9.9e-09 @ 1631 | 9.9e-09 |
| 32 | 0.0404 | 5.747e-01 | 1.548e+00 | 1.957e-01 | 2.218e+00 | 1.465e-01 | 3198 | 1.0e-08 @ 3198 | 1.0e-08 |

| N | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|
| 16, 24, 32 | 1.15 ± 0.02 | 1.13 ± 0.03 | 1.64 ± 0.02 | 0.88 ± 0.03 | 2.08 ± 0.05 |

## Serial and decomposed runs

A serial run and a four-rank scotch run (`parallel`) of the case defaults at
`N = 12` agree to a maximum relative difference over all rows of 7.5e-10.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
