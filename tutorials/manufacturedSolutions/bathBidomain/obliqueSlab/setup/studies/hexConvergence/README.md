# hexConvergence - bathBidomain/obliqueSlab

Spatial convergence of the level-set transmission solution on the
block-structured slab: orthogonal (`-shear 0`) and sheared (default shears)
hexahedra, elliptic-only and full-time.

## Ladder

| N | cells | full-time deltaT |
|---|---|---|
| 10 | 200 | Allrun default |
| 20 | 800 | Allrun default |
| 40 | 3200 | Allrun default |
| 80 | 12800 | Allrun default |

## Execution

In a copy of the case, for each `N`:

```bash
./Allrun parallel -N <n> -shear 0 -elliptic -nOuter 10 -nNonOrth 10
./Allrun parallel -N <n> -elliptic -nOuter 10 -nNonOrth 10
./Allrun parallel -N <n>
./Allrun parallel -N <n> -config ATI
```

Rows: `phiE_heart`, `phiE_bath`, `phiI`, `Vm`, `interfaceFlux`, the
`_bulk`/`_boundary` parts and the `gradVm`/`gradReconVm` boundary rows.

## Expected behaviour

Observed L2 order between successive levels (N = 10, 20, 40, 80 elliptic;
10, 20, 40 full-time; ATBI unless stated):

| run | field | orders |
|---|---|---|
| orthogonal, elliptic | `phiE_heart`, `phiI`; `interfaceFlux` | 2.02, 1.99, 2.00; 1.89, 1.88, 1.94 |
| sheared, elliptic | `phiE_heart`; `phiE_bath`; `interfaceFlux` | 1.44, 1.75, 1.91; 1.67, 1.70, 1.86; 1.18, 1.48, 1.78 |
| sheared, elliptic | `gradReconVm_boundary` | 0.60, 0.51, 0.50 |
| sheared, full-time | `Vm`; `Vm_boundary` | 1.61, 1.87; 1.68, 2.08 |
| sheared, full-time, ATI | `Vm`; `Vm_bulk` | 0.79, 0.95; 0.96, 0.59 |

On the sheared slab the exact-field wall gradient (`gradReconVm_boundary`)
converges at order one half in L2: the heart sub-mesh places the exact wall
value at the normal foot of a sheared wall cell, an O(1) error per wall cell.
With ATI the zeroGradient `Vm` wall trace lowers the `Vm` order in the bulk as
well as at the wall.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
