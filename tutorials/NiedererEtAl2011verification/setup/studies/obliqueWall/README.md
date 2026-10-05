# Oblique-wall study

This study tests the wall treatment of the monodomain solver when the
conductivity is not aligned with the walls. The Niederer et al. (2011) slab
has a conductivity tensor `diag(σl, σt, σt)` along the slab axes, so at every
wall `σn` is parallel to `n`: a zero normal gradient and a zero normal flux
are the same condition, and `sealedHeartBoundary` and `sealedWallTrace`
change nothing. The aligned slab cannot tell the treatments apart. Here the
tensor is rotated by 45° in the x–y plane, so at the x and y walls the
insulated condition `n·σ∇V = 0` needs a non-zero normal gradient, and the
treatments differ. The z walls stay aligned.

The case is the native `NiedererEtAl2011verification` case, hex mesh, with
only the patches listed below.

## What differs from `cartesianConvergence`

Each study file is `../cartesianConvergence/sweep_hex_convergence.json` with
the same grid and extra keys in `base`, and a longer `endTime` at Δx 0.1 mm.

| key (under `constant/electroProperties:monodomainSolverCoeffs.`) | value |
|---|---|
| `conductivity` | `[-1 -3 3 0 0 2 0] (0.07551194956 0.05790577195 0 0.07551194956 0 0.01760617761)`, which is `R(45°) diag(σl, σt, σt) Rᵀ` of the native `diag(0.1334177215, 0.01760617761, 0.01760617761)`, rounded to 11 digits |

Wall variants:

| variant | `sealedHeartBoundary` | `sealedWallTrace` | meaning |
|---|---|---|---|
| `0` | `false` | `zeroGradient` | the flux through the wall faces is not removed; the wall trace of `Vm` has a zero normal gradient |
| `A` | `true` | `zeroGradient` | the insulated face conductivity removes the flux through the wall faces; the wall trace keeps the zero normal gradient |
| `AB` | `true` | `conormal` | as `A`, and the wall trace is `conormalZeroFlux`, the zero-flux condition `n·σ∇V = 0` |

Time schemes, paired:

| scheme | `monodomainSolverCoeffs.timeCouplingScheme` | `system/fvSchemes:ddtSchemes.default` |
|---|---|---|
| `godunov` | `godunov` | `Euler` |
| `sbdf2` | `sbdf2` | `backward` |

`timeCouplingScheme` is read from the `monodomainSolverCoeffs` block, defaults
to `godunov` and accepts exactly these two words. `sbdf2` needs an ionic model
that supports both second-order couplings (`TNNPcompactBatched`, the case's
model, does). The diffusion step uses `fvm::ddt` with whatever scheme
`fvSchemes` names; the solver does not check it. The pairing is the study's
choice: the first-order splitting is run with the first-order ddt, and the
second-order splitting with the second-order ddt, so that neither result
mixes the orders.

The six files are `sweep_hex_oblique_<scheme>_<variant>.json`, one per
scheme and variant. Within one scheme they differ only in the two wall keys;
within one variant only in the two time-scheme keys.

## Grid and end times

Niederer's grid, as in `cartesianConvergence`: Δx = 0.5, 0.2, 0.1 mm × Δt =
0.05, 0.01, 0.005 ms, nine cases per file, the Δt values of each Δx in
order. `endTime` is 0.2 s at Δx 0.5 mm, 0.08 s at Δx 0.2 mm and 0.075 s at
Δx 0.1 mm.

The last differs from what the aligned slab needs. At Δx 0.1 mm the aligned
far corner (P8) activates at about 45.5 ms, so the 0.055 s of
`cartesianConvergence` leaves about 9.5 ms. The rotated tensor and the wall
treatment delay the last activation: at Δx 0.2 mm the last `AB` cell
activates at about 70 ms, 15.6 ms after the aligned P8 (54.4 ms). The same
delay at Δx 0.1 mm puts the last cell near 61 ms, and a cell not activated
by `endTime` cannot be compared. 0.075 s keeps about 14 ms. The run's summary
must still confirm that every cell activated.

Hex only. A tet version of the study is a later step.

## Expected outputs

Per case, as for any `niederer2011` case:

- `postProcessing/Niedererpoints/0/activationTime`: the activation times at
  P1 to P9;
- `postProcessing/Niedererlines/0/activationTime`: the activation times along
  the diagonal line;
- `<endTime>/activationTime`: the activation time of every cell, with `-1`
  where a cell did not activate;
- `<endTime>/Vm`: its wall patches carry `zeroGradient` for variants `0` and
  `A`, and `conormalZeroFlux` for `AB`.

The quantities to read are P8 against Δx and Δt per scheme and variant, the
field differences between variants (`A` − `0`, `AB` − `0`, `AB` − `A`), and
whether the `AB` − `0` gap shrinks, stays or grows as Δx goes from 0.5 to
0.1 mm.

## How to run

```bash
[omnidriver command to run]
```

The runbook for a cluster (the environment, the 18 jobs, the Slurm template,
what to copy back) is in omniD's `benchmarks/niederer2011/campaign/README.md`,
section "Oblique-wall study".
