# hexConvergence - bidomain/sphericalShell

Spatial convergence on the cubed sphere (`Allrun -mesh hex`, default): sealed
curved walls approximated by planar faces, two unequal oblique full
conductivity tensors and nearly orthogonal wall cells.

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

`N` is a multiple of 4 and `k = N/4`. The verifier writes
`postProcessing/3D_hex_<N>_cells.dat`: rows `Vm`, `u1`, `u2`, `phiE_gauge`, the
bulk, boundary and gradient rows, and the header lines `# h`, `# time` and
`# cellsPerDirection`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<sealed|stock>_hex_N<n>_s1[_T2|_dthalf|_serial]`, the names that
`../summarize.py` reads: the wall gradient (`-wallGradient`), `N`, the long
window (`-endTime 2`), half the default time step at `T = 0.2` (`-dt`), or a
serial run. From this folder:

```bash
mkdir -p results
for g in sealed stock
do
    for n in 8 12 16 24 32 40 48
    do
        rsync -a --exclude setup ../../../ results/${g}_hex_N${n}_s1/
        results/${g}_hex_N${n}_s1/Allrun parallel -mesh hex -N $n -wallGradient $g
    done
    for n in 16 24 32 48
    do
        rsync -a --exclude setup ../../../ results/${g}_hex_N${n}_s1_T2/
        results/${g}_hex_N${n}_s1_T2/Allrun parallel -mesh hex -N $n -wallGradient $g -endTime 2
    done
    for n in 32
    do
        dt=$(awk -v N=$n 'BEGIN { h = 0.00892857*(10/N)^2; k = int(0.2/h); if (k*h < 0.2) k++; printf "%.17g", 0.1/k }')
        rsync -a --exclude setup ../../../ results/${g}_hex_N${n}_s1_dthalf/
        results/${g}_hex_N${n}_s1_dthalf/Allrun parallel -mesh hex -N $n -wallGradient $g -dt $dt
    done
done
rsync -a --exclude setup ../../../ results/sealed_hex_N16_s1_serial/
results/sealed_hex_N16_s1_serial/Allrun -mesh hex -N 16
python3 ../summarize.py results
```

The time step is the case default, `endTime/ceil(endTime/(0.00892857
(10/N)^2))`, which scales as `h^2`. `-wallGradient stock` replaces the four
`sealedWallLeastSquares` entries of `system/fvSchemes` by `leastSquares`.

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`. Errors are L2 norms against the exact fields, except the
Linf columns; `h` is the verifier's `# h`. `phiE_gauge` compares `phiE` after
removing the volume-weighted mean error. The boundary rows are normalised by
the total volume, so their order is half an order above the RMS order over the
wall band. Orders are least-squares fits of `log(error)` against `log(h)` over
the three finest levels, with their standard errors.

### `sealedWallLeastSquares` (as shipped), T = 0.2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 8 | 0.1603 | 5.906e-03 | 2.464e-02 | 1.289e-01 | 3.372e-01 | 4.826e-03 | 4.288e-01 |
| 12 | 0.1074 | 3.165e-03 | 1.352e-02 | 6.803e-02 | 1.956e-01 | 2.529e-03 | 2.805e-01 |
| 16 | 0.0807 | 1.926e-03 | 9.321e-03 | 4.161e-02 | 1.231e-01 | 1.513e-03 | 1.987e-01 |
| 24 | 0.0539 | 9.021e-04 | 5.429e-03 | 2.000e-02 | 6.051e-02 | 6.742e-04 | 1.196e-01 |
| 32 | 0.0404 | 5.149e-04 | 3.534e-03 | 1.167e-02 | 3.575e-02 | 3.633e-04 | 8.241e-02 |
| 40 | 0.0324 | 3.317e-04 | 2.628e-03 | 7.622e-03 | 2.343e-02 | 2.213e-04 | 6.127e-02 |
| 48 | 0.0270 | 2.313e-04 | 2.043e-03 | 5.364e-03 | 1.655e-02 | 1.465e-04 | 4.785e-02 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 32, 40, 48 | 1.98 ± 0.00 | 1.35 ± 0.02 | 1.92 ± 0.01 | 1.90 ± 0.00 | 2.24 ± 0.01 | 1.34 ± 0.01 |

`Vm` and `phiE` are second order in L2; the wall-cell gradient is first order
(the boundary row less half an order).

### `sealedWallLeastSquares`, T = 2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 16 | 0.0807 | 5.189e-02 | 2.100e-01 | 6.789e-02 | 2.096e-01 | 3.601e-02 | 4.473e-01 |
| 24 | 0.0539 | 2.404e-02 | 1.039e-01 | 3.264e-02 | 1.033e-01 | 1.391e-02 | 2.314e-01 |
| 32 | 0.0404 | 1.383e-02 | 6.058e-02 | 1.904e-02 | 6.082e-02 | 6.977e-03 | 1.496e-01 |
| 48 | 0.0270 | 6.286e-03 | 2.751e-02 | 8.752e-03 | 2.817e-02 | 2.605e-03 | 8.218e-02 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 24, 32, 48 | 1.94 ± 0.01 | 1.92 ± 0.02 | 1.90 ± 0.01 | 1.88 ± 0.02 | 2.42 ± 0.01 | 1.49 ± 0.01 |

### Stock `leastSquares` (control), T = 0.2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 8 | 0.1603 | 5.380e-03 | 2.480e-02 | 1.131e-01 | 3.180e-01 | 4.190e-03 | 3.081e-01 |
| 12 | 0.1074 | 2.906e-03 | 1.390e-02 | 5.695e-02 | 1.592e-01 | 2.208e-03 | 2.350e-01 |
| 16 | 0.0807 | 1.889e-03 | 9.637e-03 | 3.472e-02 | 9.505e-02 | 1.475e-03 | 1.994e-01 |
| 24 | 0.0539 | 1.093e-03 | 8.651e-03 | 1.785e-02 | 4.868e-02 | 9.052e-04 | 1.626e-01 |
| 32 | 0.0404 | 7.789e-04 | 8.407e-03 | 1.163e-02 | 3.127e-02 | 6.599e-04 | 1.421e-01 |
| 40 | 0.0324 | 6.115e-04 | 7.750e-03 | 8.580e-03 | 2.465e-02 | 5.164e-04 | 1.282e-01 |
| 48 | 0.0270 | 5.059e-04 | 7.193e-03 | 6.813e-03 | 2.036e-02 | 4.204e-04 | 1.179e-01 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 32, 40, 48 | 1.07 ± 0.01 | 0.38 ± 0.01 | 1.32 ± 0.03 | 1.06 ± 0.01 | 1.11 ± 0.01 | 0.46 ± 0.00 |

`Vm` is first order: the zero-gradient wall value enters the wall-cell
gradient, whose error does not vanish with `h`.

### Stock `leastSquares` (control), T = 2

| N | h | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|---|
| 16 | 0.0807 | 4.957e-02 | 2.216e-01 | 5.728e-02 | 1.607e-01 | 3.491e-02 | 3.881e-01 |
| 24 | 0.0539 | 2.656e-02 | 1.188e-01 | 2.945e-02 | 8.473e-02 | 1.601e-02 | 2.895e-01 |
| 32 | 0.0404 | 1.809e-02 | 9.014e-02 | 1.915e-02 | 6.179e-02 | 9.738e-03 | 2.441e-01 |
| 48 | 0.0270 | 1.128e-02 | 6.059e-02 | 1.119e-02 | 4.034e-02 | 5.169e-03 | 1.961e-01 |

| N | Vm L2 | Vm Linf | phiE_gauge L2 | phiE_gauge Linf | Vm_boundary L2 | gradVm_boundary L2 |
|---|---|---|---|---|---|---|
| 24, 32, 48 | 1.23 ± 0.05 | 0.97 ± 0.01 | 1.39 ± 0.05 | 1.07 ± 0.01 | 1.63 ± 0.05 | 0.56 ± 0.01 |

### Time step

Halving the time step at `N = 32`, T = 0.2, changes the `Vm` and `phiE_gauge`
L2 errors by less than 0.1 % with either gradient: the errors above are
spatial.

### Serial and decomposed runs

A serial run and a four-rank scotch run (`parallel`) with
`sealedWallLeastSquares` at `N = 16` agree to a largest relative difference
over all rows of 6.6e-12.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
