# tetConvergence - bathBidomain/obliqueSlab

Spatial convergence of the level-set transmission solution on gmsh
tetrahedra (`Allrun -mesh tet`): the slab periodic in `tau` and `z`, the
interfaces planar and conforming.

## Ladder

| N | lc | cells |
|---|---|---|
| 10 | 0.1 | 2460 |
| 14 | 0.0714 | 5948 |
| 20 | 0.05 | 16283 |
| 28 | 0.0357 | 43272 |

Meshes at different `N` are generated independently; repeat two levels with
`-seed 2` and `-seed 3` to see the spread of the orders.

## Execution

In a copy of the case, for each `N`:

```bash
./Allrun parallel -mesh tet -N <n> -elliptic -nOuter 10 -nNonOrth 10
./Allrun parallel -mesh tet -N <n>
./Allrun parallel -mesh tet -N <n> -config ATI
```

## Expected behaviour

| run | field | orders (successive levels) |
|---|---|---|
| elliptic, ATBI | `phiE_heart`; `phiE_bath` | 2.51, 1.57, 2.48; 1.88, 1.63, 2.12 |
| elliptic, ATBI | `interfaceFlux` | 0.72, 0.91, 1.22 |
| full-time, ATBI | `Vm`; `phiE_heart` | 1.45, 1.82; 2.39, 1.64 |
| full-time, ATI | `Vm` | 0.12, 1.43 |
| any | `gradReconVm_boundary` | about one half (0.02, 0.54, 0.58) |

The tetrahedral wall cells are skewed, so the exact-field wall gradient is
O(1) per wall cell as on the sheared hexahedra. Serial and decomposed runs of
the periodic tetrahedral slab differ by about 2e-3 relative at the `tau`
seam; decompositions agree with each other to 1e-8.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
