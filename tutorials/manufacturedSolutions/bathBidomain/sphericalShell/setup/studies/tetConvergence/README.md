# tetConvergence - bathBidomain/sphericalShell

Spatial convergence on gmsh tetrahedra (`Allrun -mesh tet`): curved
interfaces with skewed wall cells.

## Ladder

| N | lc | cells |
|---|---|---|
| 7 | 0.143 | 4.5e4 |
| 10 | 0.1 | 1.3e5 |
| 14 | 0.0714 | 3.4e5 |
| 20 | 0.05 | 1.0e6 |

Meshes at different `N` are generated independently; repeat two levels with
`-seed 2` and `-seed 3` (about 3 % spread of the errors at fixed `N`).

## Execution

```bash
./Allrun parallel -mesh tet -N <n> -decomposition scotch -np <ranks> -elliptic -nOuter 10 -nNonOrth 10
./Allrun parallel -mesh tet -N <n> -decomposition scotch -np <ranks> -endTime 0.02
```

## Expected behaviour

| run | field | orders (N = 7, 10, 14, 20) |
|---|---|---|
| elliptic, ATI | `phiE_heart`; `phiE_bath` | 2.62, 2.53, 1.97; 2.55, 2.44, 2.03 |
| elliptic, ATI | `interfaceFlux`; `interfaceFluxJumpScale` | 1.47, 1.34, 0.70; 2.02, 2.08, 1.96 |
| full-time, ATI | `Vm` | about 0 (error held by the wall cells) |
| ATBI | all | diverges: the conormal wall trace closes a lagged loop with gain 1.2 - 1.4 per correction |

The 10 x 10 elliptic iteration diverges on some meshes (`-N 14 -seed 3`) with
the stock or the jump-aware gradient; the same mesh runs with one correction
per step.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
