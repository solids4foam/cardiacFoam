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

### Outer iterations and relaxation

The outer iteration and its relaxation are described in
[Outer iteration](../../../README.md#outer-iteration). These runs measure
how many outer iterations each level needs, with the case relaxation `0.9`
and without relaxation, on a doubling ladder `N = 5, 10, 20, 40` (seed 1,
four ranks, cap raised to 30000 so that only the `1e-8` tolerance stops a
run), plus `N = 7, 14, 28` with relaxation. From this folder:

```bash
mkdir -p results
for g in sealed stock
do
    for n in 5 10 20 40
    do
        rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_relaxed/
        results/${g}_tet_N${n}_relaxed/Allrun parallel -mesh tet -N $n -seed 1 -wallGradient $g -iterations 30000
        rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_unrelaxed/
        sed -i.bak '/^relaxationFactors/,/^}/d' results/${g}_tet_N${n}_unrelaxed/system/fvSolution
        results/${g}_tet_N${n}_unrelaxed/Allrun parallel -mesh tet -N $n -seed 1 -wallGradient $g -iterations 30000
    done
done
grep -H "PIMPLE: \(not \)\?converged" results/*/log.cardiacFoam
```

The `sed` line deletes the `relaxationFactors` block and leaves the
`#includeIfPresent` lines that `Allrun` options rely on.

Outer iterations to the `1e-8` tolerance:

| N | cells | sealed, relaxation 0.9 | stock, relaxation 0.9 | sealed, no relaxation | stock, no relaxation |
|---|---|---|---|---|---|
| 5 | 2302 | 280 | 314 | 102 | 89 |
| 7 | 6460 | 710 | 732 | | |
| 10 | 18177 | 1423 | 1360 | 242 | 212 |
| 14 | 46846 | 3329 | 3650 | | |
| 20 | 137124 | 8163 | 8039 | 531 | 499 |
| 28 | 368100 | 14896 | 14700 | | |
| 40 | 1075610 | 28990 | 27823 | does not converge | does not converge |

- With relaxation every level converges, and the count grows roughly with
  `N^2`. `N = 20` and finer exceed the default cap of 8000: run them with
  `-iterations`.
- Without relaxation `N = 5` to `20` converge in a few hundred iterations, to
  the same solution: the sealed `activationTime` L2 errors differ from the
  relaxed ones by 0.11 % at `N = 10` and 0.02 % at `N = 20`, and the stock
  error by 0.04 % at `N = 20`. At `N = 40` the outer residual is `1.6e-5`
  after 100 iterations and rises, slower and slower, to `1.3e-4` by about
  1000 iterations; the runs were stopped at about 1000 (sealed) and 730
  (stock) iterations.

Larger `c0` (sealed wall gradient, seed 1; `c0` edited in
`constant/electroProperties`, for example
`sed -i.bak 's/^\( *c0 *\[[^]]*\]\) *1;/\1 5;/' constant/electroProperties`).
`psi` scales as `1/c0`, so `c0 L2` compares the errors across `c0`:

| c0 | N | relaxation | outer iterations | outer residual over the last iterations | activationTime L2 | c0 L2 |
|---|---|---|---|---|---|---|
| 1 | 10 | 0.9 | 1423 | converged to `1e-8` | 1.058 | 1.058 |
| 5 | 10 | 0.9 | 200000, cap (serial) | constant at `4.0e-5` | 0.360 | 1.80 |
| 20 | 10 | 0.9 | 3000, cap (serial) | between `4.7e-4` and `5.1e-4` | 0.116 | 2.32 |
| 5 | 10 | none | stopped at about 16000 | between `1.70e-3` and `1.75e-3` | | |
| 20 | 10 | none | stopped at about 16000 | between `1.94e-3` and `2.27e-3` | | |
| 5 | 20 | none | 20000, cap | between `9.69e-5` and `9.77e-5` | 0.186 | 0.93 |

No larger-`c0` run reaches the `1e-8` tolerance, with or without relaxation,
and raising the cap does not change that: the residual has stopped
decreasing. The errors of these rows are those of unconverged iterates.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
