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
`../summarize.py` reads: the wall gradient (`-wallGradient`), `N`, the gmsh
seed, the long window (`-endTime 2`) or half the default time step at
`T = 0.2` (`-dt`). From this folder:

```bash
mkdir -p results
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
    for n in 14 40
    do
        dt=$(awk -v N=$n 'BEGIN { h = 0.00892857*(10/N)^2; k = int(0.2/h); if (k*h < 0.2) k++; printf "%.17g", 0.1/k }')
        rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_s1_dthalf/
        results/${g}_tet_N${n}_s1_dthalf/Allrun parallel -mesh tet -N $n -seed 1 -wallGradient $g -dt $dt
    done
done
python3 ../summarize.py results
```

The time step is the case default, `endTime/ceil(endTime/(0.00892857
(10/N)^2))`, which scales as `h^2`. `-wallGradient stock` replaces the two
`sealedWallLeastSquares` entries of `system/fvSchemes` by `leastSquares`.

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`. Errors are L2 norms against the exact fields, except the
`Vm` Linf column; `h` is the verifier's `# h`. At T = 0.2 the values are means
over the three seeds; the orders are least-squares fits of `log(error)`
against `log(h)` over every seed of the three finest levels, with their
standard errors. The boundary rows are normalised by the total volume, so
their order is half an order above the RMS order over the wall band.

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

Over the seeds at one `N` the `Vm` L2 error varies by at most 4.5 %
(`N = 8`), and by at most 1.4 % from `N = 20` on.

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

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
