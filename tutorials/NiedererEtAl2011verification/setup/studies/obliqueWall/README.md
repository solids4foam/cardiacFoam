# obliqueWall

The Niederer slab with its conductivity rotated 45° in the x–y plane, so the
x and y walls see a conductivity oblique to them. The aligned slab cannot tell
the wall treatments apart; this one can.

Each file is the `cartesianConvergence` grid (Δx 0.5, 0.2, 0.1 mm × Δt 0.05,
0.01, 0.005 ms) with these patches in `base`:

| patch | values |
|---|---|
| `conductivity` | `(0.07551194956 0.05790577195 0 0.07551194956 0 0.01760617761)`, R(45°) diag(σl, σt, σt) Rᵀ of the native tensor |
| wall: `0`, `A`, `AB` | `sealedHeartBoundary` false / true / true; `sealedWallTrace` zeroGradient / zeroGradient / conormal |
| scheme: `godunov`, `sbdf2` | `timeCouplingScheme` with `ddtSchemes default` Euler / backward |

`endTime` is 0.075 s at Δx 0.1 mm, above the native 0.055 s: the rotated
tensor delays the last activation by about 15 ms at Δx 0.2 mm.

Run each file, then tabulate its output:

```bash
[omnidriver command to run]
python3 setup/table_summary.py <sweep output directory>
```
