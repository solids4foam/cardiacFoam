# hexConvergence - bathBidomain/sphericalShell

Spatial convergence on the cubed sphere (`Allrun -mesh hex`, default):
curved interfaces with nearly orthogonal wall cells.

## Ladder

| N | k | cells |
|---|---|---|
| 12 | 3 | 9477 |
| 16 | 4 | 22464 |
| 24 | 6 | 75816 |
| 32 | 8 | 179712 |
| 48 | 12 | 6.1e5 |

## Execution

In a copy of the case, for each `N`:

```bash
./Allrun parallel -N <n> -decomposition scotch -np <ranks> -elliptic -nOuter 10 -nNonOrth 10
./Allrun parallel -N <n> -decomposition scotch -np <ranks>
./Allrun parallel -N <n> -decomposition scotch -np <ranks> -config ATBI
```

The elliptic recipe iterates the deferred corrections to convergence; check
that the final `phiE` residuals reach the solver tolerance, and report a level
that does not converge instead of its order.

## Expected behaviour

Observed L2 order between successive levels (N = 12, 16, 24, 32):

| run | field | orders |
|---|---|---|
| elliptic, ATI | `phiE_heart`; `phiE_bath` | 1.93, 1.90, 1.86; 1.86, 1.88, 1.88 |
| full-time, ATI | `Vm`; `Vm_bulk`; `Vm_boundary` | 0.86, 0.74, 0.66; 1.42, 1.59, 1.61; 0.55, 0.51, 0.53 |
| full-time, ATI | `gradVm_boundary`; `gradReconVm_boundary` | 0.20, 0.27, 0.35; 1.98, 2.13, 2.24 |
| full-time, ATBI | `Vm`; `phiE_heart`; `interfaceFlux` | 1.55, 1.63, 1.66; 2.01, 2.00, 2.00; 1.95, 1.96, 1.97 |
| full-time, ATBI | `interfacePartCurrent_1`, `_2`; `gradVm_boundary` | 2.57, 1.93, 2.01; 2.18, 2.03, 1.99; 1.72, 1.72, 1.66 |

With ATI the exact-field wall gradient converges at second order while the
solved one does not: the zeroGradient `Vm` wall trace, not the reconstruction,
holds `Vm` at the wall to an O(1) error per wall cell. The conormal trace
(ATBI) is consistent and converges on this mesh.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
