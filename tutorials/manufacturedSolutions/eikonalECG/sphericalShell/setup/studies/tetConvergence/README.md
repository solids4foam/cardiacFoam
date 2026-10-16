# tetConvergence - eikonalECG/sphericalShell

Spatial convergence and outer-iteration behaviour on gmsh tetrahedra
(`Allrun -mesh tet`): sealed curved walls approximated by planar faces, an
oblique full conductivity tensor and skewed wall cells.

## Ladder

| N | lc | cells (seed 1) | pinned cells |
|---|---|---|---|
| 5 | 0.2 | 2302 | 7 |
| 6 | 0.167 | 4430 | 15 |
| 7 | 0.143 | 6460 | 20 |
| 8 | 0.125 | 9815 | 37 |
| 10 | 0.1 | 18177 | 71 |
| 14 | 0.0714 | 46846 | 181 |
| 20 | 0.05 | 137124 | 505 |
| 28 | 0.0357 | 368100 | 1386 |
| 40 | 0.025 | 1075610 | 4136 |

Meshes at different `N` are generated independently, with `-seed 1`. The
levels contain the doubling ladders `5, 10, 20, 40` and `7, 14, 28`. The
verifier writes `postProcessing/3D_tet_<N>_cells.dat` (rows and header lines
as in [hexConvergence](../hexConvergence/README.md)).

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<sealed|stock>_tet_N<n>_s1`, the names that `../summarize.py` reads: the
wall gradient (`-wallGradient`), `N` and the gmsh seed. The cap is raised to
30000 so that only the `1e-8` tolerance stops a run
([Outer iterations and relaxation](#outer-iterations-and-relaxation)). From
this folder:

```bash
mkdir -p results
for g in sealed stock
do
    for n in 5 6 7 8 10 14 20 28 40
    do
        rsync -a --exclude setup ../../../ results/${g}_tet_N${n}_s1/
        results/${g}_tet_N${n}_s1/Allrun parallel -mesh tet -N $n -seed 1 -wallGradient $g -iterations 30000
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
that of the last solve. The first order table gives the pairwise order between
each level `N` and `2N`, with `r` the ratio of their `h`; the second the
least-squares fit of `log(error)` against `log(h)` over the three finest
levels, with its standard error. The rows at `N = 6` and `8` are serial runs,
the others four-rank runs.

With the case defaults (advection-diffusion formulation, relaxation `0.9`)
every level of both wall gradients converges to the `1e-8` outer tolerance.
With `sealedWallLeastSquares` the `activationTime` L2 order rises along both
doubling ladders and reaches 2.00 (Linf 1.99) at `N = 20` to `40`; the coarser
pairs (1.17 to 1.83) are not yet asymptotic, and fits over short ladders
inside that range understate the order. The stock `leastSquares` gradient
falls back to 1.25 (Linf 1.14) at `N = 20` to `40` and 1.30 at `N = 14` to
`28`. At `N = 40` the stock L2 error is still marginally the smaller one
(7.6e-2 against 8.3e-2), while the sealed Linf error is the smaller one
(0.20 against 0.30).

### sealed, tet, advectionDiffusion, converged

| N | h | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 | outer its | min residual | last residual |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.1164 | 2.356e+00 | 6.192e+00 | 1.451e+00 | 4.604e+00 | 3.292e+00 | 280 | 9.3e-09 @ 280 | 9.3e-09 |
| 6 | 0.0937 | 1.918e+00 | 4.686e+00 | 1.016e+00 | 4.140e+00 | 3.107e+00 | 617 | 9.6e-09 @ 617 | 9.6e-09 |
| 7 | 0.0827 | 1.615e+00 | 4.283e+00 | 8.022e-01 | 3.873e+00 | 3.098e+00 | 710 | 9.7e-09 @ 710 | 9.7e-09 |
| 8 | 0.0719 | 1.289e+00 | 3.354e+00 | 5.952e-01 | 3.594e+00 | 2.927e+00 | 905 | 9.8e-09 @ 905 | 9.8e-09 |
| 10 | 0.0586 | 1.058e+00 | 2.594e+00 | 4.313e-01 | 3.088e+00 | 2.610e+00 | 1423 | 9.8e-09 @ 1423 | 9.8e-09 |
| 14 | 0.0428 | 6.616e-01 | 1.606e+00 | 2.348e-01 | 2.428e+00 | 2.034e+00 | 3329 | 1.0e-08 @ 3329 | 1.0e-08 |
| 20 | 0.0299 | 3.270e-01 | 7.919e-01 | 9.886e-02 | 1.709e+00 | 1.423e+00 | 8163 | 1.0e-08 @ 8163 | 1.0e-08 |
| 28 | 0.0215 | 1.882e-01 | 4.165e-01 | 4.858e-02 | 1.166e+00 | 9.648e-01 | 14896 | 1.0e-08 @ 14896 | 1.0e-08 |
| 40 | 0.0150 | 8.300e-02 | 2.014e-01 | 1.880e-02 | 7.508e-01 | 6.151e-01 | 28990 | 1.0e-08 @ 28990 | 1.0e-08 |

| N | r | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|---|
| 5 to 10 | 1.987 | 1.17 | 1.27 | 1.77 | 0.58 | 0.34 |
| 7 to 14 | 1.933 | 1.35 | 1.49 | 1.86 | 0.71 | 0.64 |
| 10 to 20 | 1.960 | 1.75 | 1.76 | 2.19 | 0.88 | 0.90 |
| 14 to 28 | 1.988 | 1.83 | 1.97 | 2.29 | 1.07 | 1.09 |
| 20 to 40 | 1.987 | 2.00 | 1.99 | 2.42 | 1.20 | 1.22 |

| N | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|
| 20, 28, 40 | 2.00 ± 0.18 | 1.99 ± 0.02 | 2.42 ± 0.14 | 1.20 ± 0.02 | 1.22 ± 0.02 |

### stock, tet, advectionDiffusion, converged

| N | h | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 | outer its | min residual | last residual |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.1164 | 1.879e+00 | 5.092e+00 | 1.156e+00 | 8.978e+00 | 3.313e+00 | 314 | 9.5e-09 @ 314 | 9.5e-09 |
| 6 | 0.0937 | 1.364e+00 | 3.757e+00 | 7.177e-01 | 7.064e+00 | 2.768e+00 | 611 | 9.6e-09 @ 611 | 9.6e-09 |
| 7 | 0.0827 | 8.920e-01 | 2.668e+00 | 4.636e-01 | 6.014e+00 | 2.434e+00 | 732 | 9.8e-09 @ 732 | 9.8e-09 |
| 8 | 0.0719 | 8.455e-01 | 2.367e+00 | 4.070e-01 | 5.298e+00 | 2.107e+00 | 1017 | 9.9e-09 @ 1017 | 9.9e-09 |
| 10 | 0.0586 | 6.440e-01 | 1.907e+00 | 2.713e-01 | 4.311e+00 | 1.779e+00 | 1360 | 9.8e-09 @ 1360 | 9.8e-09 |
| 14 | 0.0428 | 3.321e-01 | 1.001e+00 | 1.336e-01 | 3.188e+00 | 1.332e+00 | 3650 | 1.0e-08 @ 3650 | 1.0e-08 |
| 20 | 0.0299 | 1.802e-01 | 6.473e-01 | 6.442e-02 | 2.392e+00 | 1.022e+00 | 8039 | 1.0e-08 @ 8039 | 1.0e-08 |
| 28 | 0.0215 | 1.356e-01 | 5.003e-01 | 3.798e-02 | 1.892e+00 | 8.293e-01 | 14700 | 1.0e-08 @ 14700 | 1.0e-08 |
| 40 | 0.0150 | 7.621e-02 | 2.956e-01 | 2.016e-02 | 1.508e+00 | 6.737e-01 | 27823 | 1.0e-08 @ 27823 | 1.0e-08 |

| N | r | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|---|
| 5 to 10 | 1.987 | 1.56 | 1.43 | 2.11 | 1.07 | 0.91 |
| 7 to 14 | 1.933 | 1.50 | 1.49 | 1.89 | 0.96 | 0.91 |
| 10 to 20 | 1.960 | 1.89 | 1.61 | 2.14 | 0.88 | 0.82 |
| 14 to 28 | 1.988 | 1.30 | 1.01 | 1.83 | 0.76 | 0.69 |
| 20 to 40 | 1.987 | 1.25 | 1.14 | 1.69 | 0.67 | 0.61 |

| N | activationTime L2 | activationTime Linf | activationTime_boundary L2 | gradActivationTime_boundary L2 | gradReconActivationTime_boundary L2 |
|---|---|---|---|---|---|
| 20, 28, 40 | 1.26 ± 0.22 | 1.15 ± 0.20 | 1.69 ± 0.05 | 0.67 ± 0.02 | 0.61 ± 0.02 |

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
