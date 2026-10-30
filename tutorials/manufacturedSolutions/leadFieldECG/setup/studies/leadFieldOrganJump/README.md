# leadFieldOrganJump - leadFieldECG

With torso layers of different conductivity and a current across them
(`-problem layered`), at what order do the lead field and the direct forward
solve converge to the exact potential? And is the ECG of the insulated case
(`-alpha 0`) independent of the layer conductivities?

The problem is described in the case [README](../../../README.md): torso
layers `sigma_1 = 0.1` (next to the heart) and `sigma_2 = 0.02` (beyond
`x = -0.5` and `x = 1.5`), a jump of 5, and the surface currents `-alpha` on
`xMin` and `alpha` on `xMax` with `alpha = 0.01`. The torso faces between the
layers take the distance-weighted harmonic mean of the two conductivities,
and the heart-torso faces the flux-continuous face conductivity of the
`distanceWeightedHarmonic` interface option, in the forward operator and
therefore in the lead fields.

The rows are the static rows of the verifier (the exact `Vm` at the heart
cell centres) and the reciprocity gap, with the electrode groups of
[leadFieldConvergence](../leadFieldConvergence/README.md): `torso`
(`torsoLeft`, `torsoLeftFar`, `torsoRight`) and `near` (`nearLeft`,
`nearRight`). Each entry is the largest Linf over the samples among the
electrodes of the group.

The problem is linear in its data, so the run with `alpha = 0.01` minus the
run with `-alpha 0` on the same mesh is the response to the current alone.
Its error, the **alpha part**, is reported for two groups:

- `left` (`torsoLeft`, `nearLeft`, in the inner left layer): between the
  pinned outer layer and these electrodes the current crosses the layer
  interface `x = -0.5` and no heart interface;
- `right` (`torsoRight`, `nearRight`): the current also crosses both
  heart-torso interfaces and the layer interface `x = 1.5`.

`torsoLeftFar` lies in the pinned layer and is in neither group.

## Ladder

| mesh | N | cells |
|---|---|---|
| hex | 10, 20, 40 | 3000, 24000, 192000 |
| tet (gmsh 4.15.2, seed 1, `layered_box.geo.template`) | 10, 14, 20, 28 | 15443, 39418, 112979, 302837 |

All runs are serial, monodomain with `corrected` snGrad, and use
`-endTime 0.02` (two samples, `dt = 0.01`). The fixed points stop at the
`leadFieldFixedPoint` tolerance of `system/fvSolution`. On hexahedra the
lead fields need 5 or 6 passes and each direct solve 2. On tetrahedra the
lead fields need 40 to 64 passes and each direct solve 39 to 77.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<configuration>_<mesh>_N<n>`, the names that `../summarize.py` reads and
pairs: `layered` with `layeredAlpha0` for the alpha part, `layeredAlpha0`
with `slab` for the invariance row. From this folder:

```bash
mkdir -p results
run()
{
    name=$1; shift
    rsync -a --exclude setup ../../../ results/$name/
    bash results/$name/Allrun "$@" -endTime 0.02
}
for n in 10 20 40
do
    run layered_hex_N$n -problem layered -mesh hex -N $n
    run layeredAlpha0_hex_N$n -problem layered -alpha 0 -mesh hex -N $n
    run slab_hex_N$n -problem slab -mesh hex -N $n
done
for n in 10 14 20 28
do
    run layered_tet_N$n -problem layered -mesh tet -N $n
    run layeredAlpha0_tet_N$n -problem layered -alpha 0 -mesh tet -N $n
    run slab_tet_N$n -problem slab -mesh tet -N $n
done
python3 ../summarize.py results
```

## Measured behaviour

Errors against the exact potential, as the largest Linf over the samples
among the electrodes of the group; the order is between consecutive levels,
with `h = (cells/3)^(-1/3)`, which is `1/N` on hexahedra. Gaps are relative to
the largest exact potential (about 1.77 with `alpha = 0.01`).

### Hexahedra

Reciprocity gap `leadField-direct` (on `Vm_h`) and static, both groups, for
`layered` and `layeredAlpha0`: at most `4.6e-12` (N = 10 to 40); the lead
field and the direct solve agree to every printed digit below.

Static errors at the cell centroid (`torso` and `near` give the same values,
set by the right-slab electrodes):

| N | LF torso | order | direct torso | order | `-alpha 0` torso | order |
|---|---|---|---|---|---|---|
| 10 | 1.7585e-02 |  | 1.7585e-02 |  | 1.7585e-02 |  |
| 20 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 |
| 40 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 |

Alpha part:

| N | LF left | LF right | direct left | direct right |
|---|---|---|---|---|
| 10 | 7.7383e-14 | 5.1292e-14 | 7.0360e-13 | 7.0588e-13 |
| 20 | 6.3405e-13 | 3.2152e-13 | 4.7679e-13 | 4.7740e-13 |
| 40 | 2.7947e-12 | 3.5263e-12 | 7.1154e-12 | 7.1378e-12 |

Static errors at the electrode point (`torso`): 1.4235e-02, 1.3553e-02,
4.4000e-03 (orders 0.07, 1.62). With `-alpha 0` the potential is constant
in each torso layer and the point errors equal the centroid errors. With
`alpha != 0` the potential has a slope of up to 0.5 in the outer layers, so
the point error contains the localisation `phi(x_e) - phi(x_c)` of the
electrode in its cell, which is `O(h)` and depends on where the electrode
falls in its cell.

- **The interfaces are exact on the hexahedra.** The alpha part is round-off
  (`<= 7.2e-12`) at both groups, in the lead field and the direct solve:
  - `left`: the harmonic face conductivity of the torso layer interface
    `x = -0.5` gives the exact flux of a piecewise linear potential across a
    jump of 5;
  - `right`: the current also crosses both heart-torso walls, whose heart
    side carries `Gi + Ge`, the conductivity of the total current through a
    wall where the sealed `Vm` carries no flux.
- **The total error is the slab's.** The layered static errors equal the
  `-alpha 0` column, second order, so the layers and the current add no
  discretisation error on these meshes.

### Tetrahedra

| N | h | gap dyn | gap static | LF torso | order | direct torso | order | LF near | direct near |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.05792 | 2.029e-10 | 2.032e-10 | 6.3479e-03 |  | 6.3479e-03 |  | 6.7079e-03 | 6.7079e-03 |
| 14 | 0.04238 | 1.460e-10 | 1.536e-10 | 3.4584e-03 | 1.94 | 3.4584e-03 | 1.94 | 3.4239e-03 | 3.4239e-03 |
| 20 | 0.02983 | 2.977e-11 | 3.153e-11 | 1.6333e-03 | 2.14 | 1.6333e-03 | 2.14 | 1.6639e-03 | 1.6639e-03 |
| 28 | 0.02148 | 4.488e-11 | 4.432e-11 | 9.0280e-04 | 1.80 | 9.0280e-04 | 1.80 | 9.2805e-04 | 9.2805e-04 |

(`gap`: the `torso` reciprocity gap on `Vm_h` and on the exact `Vm`.)

Static errors at the electrode point, lead field and direct: `torso`
2.1574e-02, 1.5392e-02, 2.7442e-03, 4.7860e-03 (orders 1.08, 4.91, -1.69);
`near` 6.0150e-03, 1.9968e-03, 6.6563e-04, 5.3771e-04 (orders 3.53, 3.13,
0.65). With `-alpha 0` they equal the centroid errors.

Alpha part (the lead field and the direct solve give the same values):

| N | left | order | right | order |
|---|---|---|---|---|
| 10 | 4.4046e-04 |  | 4.4610e-04 |  |
| 14 | 2.2056e-04 | 2.21 | 3.1302e-04 | 1.13 |
| 20 | 4.2193e-05 | 4.71 | 9.9369e-05 | 3.27 |
| 28 | 2.9254e-05 | 1.11 | 9.7180e-05 | 0.07 |

- **Lead field and direct solve.** The reciprocity gap is at the level of the
  fixed-point tolerances (`<= 2.6e-10`), and the lead field equals the
  direct solve in every row, alpha part included.
- **Total error.** The centroid errors converge at about second order with
  scatter, 1.80 to 2.14 in the torso group.
- **Point errors.** The point-evaluated electrode errors converge with
  erratic level-to-level orders, in both groups.
- **Alpha part.** It falls from 4.4e-4 to 2.9e-5 (`left`) and to 9.7e-5
  (`right`), with erratic level-to-level orders: `left` 2.21, 4.71, 1.11;
  `right` 1.13, 3.27, 0.07. No order is measured, and the cause of the
  scatter is not established.

### Invariance: `-alpha 0` against the slab

ECG `torsoRight - torsoLeft`, largest difference over the samples between
`layeredAlpha0` and `slab` at the same `N`, from the channel on `Vm_h` and
from the static channel, and each channel's error against the exact
`2 k S`:

| mesh | N | difference (Vm_h) | difference (static) | error layered | error slab |
|---|---|---|---|---|---|
| hex | 10 | 1.503e-13 | 1.510e-13 | 1.7645e-02 | 1.7645e-02 |
| hex | 20 | 2.022e-12 | 2.005e-12 | 4.4945e-03 | 4.4945e-03 |
| hex | 40 | 2.053e-12 | 2.044e-12 | 1.2006e-03 | 1.2006e-03 |
| tet | 10 | 9.078e-05 | 6.849e-05 | 6.7147e-03 | 6.6239e-03 |
| tet | 14 | 1.002e-04 | 7.394e-05 | 3.6883e-03 | 3.5881e-03 |
| tet | 20 | 1.754e-05 | 2.175e-05 | 1.9254e-03 | 1.9396e-03 |
| tet | 28 | 2.723e-05 | 2.050e-05 | 1.0550e-03 | 1.0278e-03 |

- **Hexahedra.** The two cases have the same cells (five blocks against
  three, the same `h`), and the ECG does not move with the layer
  conductivities: the difference is round-off.
- **Tetrahedra.** The layered mesh has two more internal planes, so the two
  cases are different meshes. Both ECGs converge to `2 k S`, and their
  difference, 1e-4 to 2e-5, is 1 to 3% of either error.

### Linear interpolation of the layer conductivity

Torso faces between different conductivities take the distance-weighted
harmonic mean. On the hexahedra it is exact for the piecewise linear
potential (the `left` rows). With the linear mean `(sigma_1 + sigma_2)/2`
instead, the face between the layers would carry the flux of a resistance
`h/(sigma_1 + sigma_2) * 2` instead of `h/2 (1/sigma_1 + 1/sigma_2)`, an
offset of `alpha h ((1/sigma_1 + 1/sigma_2)/2 - 2/(sigma_1 + sigma_2))` per
interface (`13.3 alpha h` here): first order at the jump.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
