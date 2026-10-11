# tetConvergence - eikonalECG/sphericalShell

Spatial convergence and outer-iteration behaviour on gmsh tetrahedra
(`Allrun -mesh tet`): sealed curved walls approximated by planar faces, an
oblique full conductivity tensor and skewed wall cells.

## Ladder

| N | lc | cells (seed 1) |
|---|---|---|
| 6 | 0.167 | 4430 |
| 8 | 0.125 | 9815 |
| 10 | 0.1 | 18177 |
| 14 | 0.0714 | 46846 |

Meshes at different `N` are generated independently, with `-seed 1`. The
verifier writes `postProcessing/3D_tet_<N>_cells.dat` (rows and header lines
as in [hexConvergence](../hexConvergence/README.md)).

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<sealed|stock>_tet_N<n>_s1`, the names that `../summarize.py` reads: the
wall gradient (`-wallGradient`), `N` and the gmsh seed. The runs are serial.
From this folder:

```bash
mkdir -p results
for g in sealed stock
do
    for n in 6 8 10 14
    do
        rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_s1/
        results/${g}_tet_N${n}_s1/Allrun -mesh tet -N $n -seed 1 -wallGradient $g
    done
done
python3 ../summarize.py results
```

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`. Errors are L2 norms against the exact fields, except the
`Linf` column; `h` is the verifier's `# h`. `outer its` is the number of outer
iterations, `min residual` the smallest initial residual of the
`activationTime` solves with the iteration where it occurs, `last residual`
that of the last solve. Orders are least-squares fits of `log(error)` against
`log(h)` over the three finest levels, with their standard errors. Both wall
gradients converge to the `1e-8` outer tolerance with the case defaults
(advection-diffusion formulation).

### sealed, tet, advectionDiffusion, converged

| N | h | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 | outer its | min residual | last residual |
|---|---|---|---|---|---|---|---|---|---|
| 6 | 0.0937 | 1.918e+00 | 4.686e+00 | 1.016e+00 | 4.140e+00 | 3.107e+00 | 617 | 9.6e-09 @ 617 | 9.6e-09 |
| 8 | 0.0719 | 1.289e+00 | 3.354e+00 | 5.952e-01 | 3.594e+00 | 2.927e+00 | 905 | 9.8e-09 @ 905 | 9.8e-09 |
| 10 | 0.0586 | 1.059e+00 | 2.594e+00 | 4.313e-01 | 3.088e+00 | 2.610e+00 | 1422 | 9.9e-09 @ 1422 | 9.9e-09 |
| 14 | 0.0428 | 6.617e-01 | 1.606e+00 | 2.348e-01 | 2.428e+00 | 2.034e+00 | 3329 | 9.9e-09 @ 3329 | 9.9e-09 |

| N | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|
| 8, 10, 14 | 1.30 ± 0.14 | 1.42 ± 0.07 | 1.80 ± 0.10 | 0.75 ± 0.01 | 0.71 ± 0.06 |

### stock, tet, advectionDiffusion, converged

| N | h | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 | outer its | min residual | last residual |
|---|---|---|---|---|---|---|---|---|---|
| 6 | 0.0937 | 1.364e+00 | 3.757e+00 | 7.177e-01 | 7.064e+00 | 2.768e+00 | 611 | 9.6e-09 @ 611 | 9.6e-09 |
| 8 | 0.0719 | 8.455e-01 | 2.367e+00 | 4.070e-01 | 5.298e+00 | 2.107e+00 | 1017 | 9.9e-09 @ 1017 | 9.9e-09 |
| 10 | 0.0586 | 6.440e-01 | 1.907e+00 | 2.713e-01 | 4.311e+00 | 1.779e+00 | 1360 | 9.9e-09 @ 1360 | 9.9e-09 |
| 14 | 0.0428 | 3.320e-01 | 1.001e+00 | 1.336e-01 | 3.188e+00 | 1.332e+00 | 3648 | 1.0e-08 @ 3648 | 1.0e-08 |

| N | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|
| 8, 10, 14 | 1.82 ± 0.21 | 1.69 ± 0.27 | 2.15 ± 0.07 | 0.97 ± 0.01 | 0.89 ± 0.03 |

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
