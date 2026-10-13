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
`<sealed|stock>_hex_N<n>_s1[_ex][_exact|_cap|_parallel]`, the names that
`../summarize.py` reads: the wall gradient (`-wallGradient`), `N`, the
explicit formulation (`_ex`, `-formulation explicit`), one iteration from the
exact field (`-exactStart -iterations 1`), an iteration cap
(`-iterations 1500`), or a parallel run. The other runs are serial. From this
folder:

```bash
mkdir -p results
for n in 8 12 16 24 32 48
do
    rsync -a --exclude setup ../../../ results/sealed_hex_N${n}_s1/
    results/sealed_hex_N${n}_s1/Allrun -mesh hex -N $n
    rsync -a --exclude setup ../../../ results/sealed_hex_N${n}_s1_exact/
    results/sealed_hex_N${n}_s1_exact/Allrun -mesh hex -N $n -exactStart -iterations 1
done
for n in 8 12 16 24 32
do
    rsync -a --exclude setup ../../../ results/stock_hex_N${n}_s1/
    results/stock_hex_N${n}_s1/Allrun -mesh hex -N $n -wallGradient stock
done
for n in 12 16
do
    rsync -a --exclude setup ../../../ results/sealed_hex_N${n}_s1_ex_cap/
    results/sealed_hex_N${n}_s1_ex_cap/Allrun -mesh hex -N $n -formulation explicit -iterations 1500
done
rsync -a --exclude setup ../../../ results/sealed_hex_N12_s1_parallel/
results/sealed_hex_N12_s1_parallel/Allrun parallel -mesh hex -N 12
python3 ../summarize.py results
```

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`. Errors are L2 norms against the exact fields, except the
`Linf` column; `h` is the verifier's `# h`. `outer its` is the number of outer
iterations, `min residual` the smallest initial residual of the
`activationTime` solves with the iteration where it occurs, `last residual`
that of the last solve. Orders are least-squares fits of `log(error)` against
`log(h)` over the three finest levels, with their standard errors.
`activationTime_boundary` is normalised by the total volume, so its order is
half an order above the RMS order over the wall band.

With the case defaults (advection-diffusion formulation) both wall gradients
converge to the `1e-8` outer tolerance; the number of outer iterations grows
roughly as `N^2` (154 at `N = 8`, 3069 at `N = 32`), and the `N = 48` run
stops at the 8000-iteration cap with the outer residual at 1.3e-8 (falling by
a factor of about 2.8 over its last 1000 iterations). With
`sealedWallLeastSquares` the `activationTime` L2 error falls at a fitted order
of 1.5 over `N = 24`-`48`, rising from 1.2 over `N = 16`-`32`. One iteration
from the exact field (`-exactStart -iterations 1`) gives an error that falls
at order 2.7 and is about a hundredth of the converged error (0.7 % to 1.6 %
over `N = 24`-`48`): most of the converged error builds up over the outer
iteration, not in one linearised solve. With the explicit formulation and
`sealedWallLeastSquares` the outer iteration does not converge: the outer
residual passes a minimum and the activation time drifts away from the pinned
cells.

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

### Serial and decomposed runs

A serial run and a four-rank scotch run (`parallel`) of the case defaults at
`N = 12` agree to a largest relative difference over all rows of 2.1e-7; both
stop at the `1e-8` outer tolerance after 303 outer iterations. Repeated
parallel runs differ from each other by up to about 1e-5 relative, within
what that tolerance allows.

### Outer iterations and relaxation

The outer iteration and its relaxation are described in
[Outer iteration](../../../README.md#outer-iteration). These runs measure
how many outer iterations each level needs, with the case relaxation `0.9`
and without relaxation, on the doubling ladders `N = 12, 24, 48` and
`N = 20, 40, 80` (four ranks, cap raised to 30000 so that only the `1e-8`
tolerance stops a run). From this folder:

```bash
mkdir -p results
for g in sealed stock
do
    for n in 12 20 24 40 48 80
    do
        rsync -a --exclude setup ../../../ results/${g}_hex_N${n}_relaxed/
        results/${g}_hex_N${n}_relaxed/Allrun parallel -mesh hex -N $n -wallGradient $g -iterations 30000
        rsync -a --exclude setup ../../../ results/${g}_hex_N${n}_unrelaxed/
        sed -i.bak '/^relaxationFactors/,/^}/d' results/${g}_hex_N${n}_unrelaxed/system/fvSolution
        results/${g}_hex_N${n}_unrelaxed/Allrun parallel -mesh hex -N $n -wallGradient $g -iterations 30000
    done
done
grep -H "PIMPLE: \(not \)\?converged" results/*/log.cardiacFoam
```

The `sed` line deletes the `relaxationFactors` block and leaves the
`#includeIfPresent` lines that `Allrun` options rely on.

Outer iterations to the `1e-8` tolerance:

| N | cells | sealed, relaxation 0.9 | stock, relaxation 0.9 | sealed, no relaxation | stock, no relaxation |
|---|---|---|---|---|---|
| 12 | 2916 | 303 | 338 | 110 | 106 |
| 20 | 13500 | 929 | 1038 | 208 | 214 |
| 24 | 23328 | 1450 | 1631 | 258 | 271 |
| 40 | 108000 | 5437 | 5236 | does not converge | does not converge |
| 48 | 186624 | 8253 | 7633 | | |
| 80 | 864000 | 21242 | 19901 | | |

- With relaxation every level converges, and the count grows roughly with
  `N^2`. The sealed `N = 48` run needs 8253 iterations, just above the
  default cap of 8000; with the cap raised its `activationTime` L2 error is
  the same to four digits as that of the capped run in the table above.
  `N = 48` and finer need `-iterations`.
- Without relaxation `N = 12` to `24` converge in a few hundred iterations,
  to the same `activationTime` L2 errors to four digits. At `N = 40` the
  outer residual falls to about `2e-2` within five iterations and then
  oscillates between `4.3e-3` and `4.8e-3`; the runs were stopped at about
  4300 (sealed) and 5000 (stock) iterations. The unrelaxed `N = 48` and
  `N = 80` levels were not run.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
