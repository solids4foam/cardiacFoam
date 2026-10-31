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
[leadFieldConvergence](../leadFieldConvergence/README.md): the body-surface
electrodes `surface` (`surface1` to `surface8`), which are the ECG and carry
the convergence result, and the interior probes `torso` (`torsoLeft`,
`torsoLeftFar`, `torsoRight`) and `near` (`nearLeft`, `nearRight`), which
show where an error comes from. Each electrode entry is the largest Linf
over the samples among the electrodes of the group; the field row
`phiStatic_surface` is the RMS and the maximum error of the direct static
potential over the cells of the outer boundary. The surface electrodes lie
in both torso layers: `surface5` and `surface7` in the inner layers, the
other six in the outer layers.

The problem is linear in its data, so the run with `alpha = 0.01` minus the
run with `-alpha 0` on the same mesh is the response to the current alone.
Its error, the **alpha part**, is reported for two groups of probes:

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
| tet (gmsh 4.15.2, seed 1, `layered_box.geo.template`), `layered` and `layeredAlpha0` | 10, 14, 20, 28, 40 | 15443, 39418, 112979, 302837, 875170 |
| tet (gmsh 4.15.2, seed 1), `slab` | 10, 14, 20, 28 | 14792, 39131, 111482, 299957 |

The tetrahedral `N = 40` runs use two ranks (the `parallel` argument of
`Allrun`); every other run is serial. All runs are monodomain with `corrected` snGrad, and use
`-endTime 0.02` (two samples, `dt = 0.01`). Each run builds 13 lead fields,
one per electrode. The fixed points stop at the `leadFieldFixedPoint`
tolerance of `system/fvSolution`. On hexahedra the lead fields need 5 or 6
passes and each direct solve 2. On tetrahedra the lead fields need 38 to 64
passes and each direct solve 39 to 77.

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
    rsync -a --exclude setup --exclude regression --exclude sphericalShell \
        ../../../ results/$name/
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
run layered_tet_N40 parallel -problem layered -mesh tet -N 40
run layeredAlpha0_tet_N40 parallel -problem layered -alpha 0 -mesh tet -N 40
python3 ../summarize.py results
```

## Measured behaviour

Errors against the exact potential, as the largest Linf over the samples
among the electrodes of the group, and the error of the direct static field
over the surface cells; the order is between consecutive levels, with
`h = (cells/3)^(-1/3)`, which is `1/N` on hexahedra. Gaps are relative to
the largest exact potential (about 1.77 with `alpha = 0.01`). The lead field
and the direct solve give the same electrode errors to every printed digit,
so the probe and point tables give the lead field only.

### Hexahedra

Reciprocity gap `leadField-direct` (on `Vm_h`) and static, every group, for
`layered` and `layeredAlpha0`: at most 4.9e-12 (N = 10 to 40).

Body surface, static errors at the cell centroid, with the `-alpha 0`
lead field:

| N | LF surface | order | direct surface | order | `phiStatic_surface` RMS | order | `phiStatic_surface` max | order | `-alpha 0` surface | order |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 1.7585e-02 |  | 1.7585e-02 |  | 1.1470e-02 |  | 1.7585e-02 |  | 1.7585e-02 |  |
| 20 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 | 2.8772e-03 | 2.00 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 |
| 40 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 | 7.2023e-04 | 2.00 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 |

The `torso` and `near` probes give the surface values to every printed
digit, set by the right-slab electrodes.

Alpha part:

| N | LF left | LF right | direct left | direct right |
|---|---|---|---|---|
| 10 | 7.7383e-14 | 5.1292e-14 | 7.0360e-13 | 7.0588e-13 |
| 20 | 6.3405e-13 | 3.2152e-13 | 4.7679e-13 | 4.7740e-13 |
| 40 | 2.7947e-12 | 3.5263e-12 | 7.1154e-12 | 7.1378e-12 |

Static errors at the electrode point:

| N | LF surface | order | LF torso | order | LF near | order |
|---|---|---|---|---|---|---|
| 10 | 4.1085e-02 |  | 1.4235e-02 |  | 1.4455e-02 |  |
| 20 | 1.5403e-02 | 1.42 | 1.3553e-02 | 0.07 | 3.7729e-03 | 1.94 |
| 40 | 5.8512e-03 | 1.40 | 4.4000e-03 | 1.62 | 1.7212e-03 | 1.13 |

With `-alpha 0` the potential is constant in each torso layer and the point
errors equal the centroid errors. With `alpha != 0` the potential has a
slope of up to 0.5 in the outer layers, so the point error contains the
localisation `phi(x_e) - phi(x_c)` of the electrode in its cell, which is
`O(h)` and depends on where the electrode falls in its cell.

- **The total error is the slab's.** The body-surface ECG and the error over
  the surface cells converge at second order and equal the `-alpha 0`
  values, so the layers and the current add no discretisation error on
  these meshes.
- **The interfaces are exact on the hexahedra.** The alpha part is round-off
  (`<= 7.2e-12`) at both probe groups, in the lead field and the direct
  solve:
  - `left`: the harmonic face conductivity of the torso layer interface
    `x = -0.5` gives the exact flux of a piecewise linear potential across a
    jump of 5;
  - `right`: the current also crosses both heart-torso walls, whose heart
    side carries `Gi + Ge`, the conductivity of the total current through a
    wall where the sealed `Vm` carries no flux.
- **Point errors.** On the body surface the point errors converge steadily
  at an order between one and two, between the first order of the
  localisation and the second order of the solve; the single probes give
  erratic orders.

### Tetrahedra

Body surface, static errors at the cell centroid, and the reciprocity gap of
the surface group on `Vm_h` (`gap dyn`) and on the exact `Vm`
(`gap static`):

| N | h | gap dyn | gap static | LF surface | order | direct surface | order | `phiStatic_surface` RMS | order | `phiStatic_surface` max | order |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.05792 | 1.871e-10 | 1.862e-10 | 6.4057e-03 |  | 6.4057e-03 |  | 4.0486e-03 |  | 1.0375e-02 |  |
| 14 | 0.04238 | 1.574e-10 | 1.652e-10 | 3.2807e-03 | 2.14 | 3.2807e-03 | 2.14 | 2.1773e-03 | 1.99 | 6.4399e-03 | 1.53 |
| 20 | 0.02983 | 2.578e-11 | 2.632e-11 | 1.7302e-03 | 1.82 | 1.7302e-03 | 1.82 | 1.1382e-03 | 1.85 | 4.5193e-03 | 1.01 |
| 28 | 0.02148 | 4.978e-11 | 4.917e-11 | 9.5604e-04 | 1.80 | 9.5604e-04 | 1.80 | 6.2835e-04 | 1.81 | 2.5761e-03 | 1.71 |
| 40 | 0.01508 | 2.385e-11 | 2.388e-11 | 4.7978e-04 | 1.95 | 4.7978e-04 | 1.95 | 3.1108e-04 | 1.99 | 1.8483e-03 | 0.94 |

Body surface with `-alpha 0`:

| N | LF surface | order | `phiStatic_surface` RMS | order | `phiStatic_surface` max | order |
|---|---|---|---|---|---|---|
| 10 | 5.8769e-03 |  | 3.7962e-03 |  | 5.8809e-03 |  |
| 14 | 3.1477e-03 | 2.00 | 2.0408e-03 | 1.99 | 3.2037e-03 | 1.94 |
| 20 | 1.5821e-03 | 1.96 | 1.0324e-03 | 1.94 | 1.5850e-03 | 2.01 |
| 28 | 8.3668e-04 | 1.94 | 5.4473e-04 | 1.95 | 8.4103e-04 | 1.93 |
| 40 | 4.0991e-04 | 2.02 | 2.6720e-04 | 2.01 | 4.2167e-04 | 1.95 |

Probes, `torso` reciprocity gap and lead-field errors:

| N | gap dyn | gap static | LF torso | order | LF near | order |
|---|---|---|---|---|---|---|
| 10 | 1.799e-10 | 1.801e-10 | 6.3479e-03 |  | 6.7079e-03 |  |
| 14 | 1.294e-10 | 1.361e-10 | 3.4584e-03 | 1.94 | 3.4239e-03 | 2.15 |
| 20 | 2.640e-11 | 2.796e-11 | 1.6333e-03 | 2.14 | 1.6639e-03 | 2.06 |
| 28 | 3.978e-11 | 3.929e-11 | 9.0280e-04 | 1.80 | 9.2805e-04 | 1.78 |
| 40 | 1.148e-11 | 2.468e-11 | 4.6697e-04 | 1.86 | 4.0645e-04 | 2.33 |

Static errors at the electrode point (with `-alpha 0` they equal the
centroid errors):

| N | LF surface | order | LF torso | order | LF near | order |
|---|---|---|---|---|---|---|
| 10 | 2.3554e-02 |  | 2.1574e-02 |  | 6.0150e-03 |  |
| 14 | 1.1292e-02 | 2.35 | 1.5392e-02 | 1.08 | 1.9968e-03 | 3.53 |
| 20 | 7.7774e-03 | 1.06 | 2.7442e-03 | 4.91 | 6.6563e-04 | 3.13 |
| 28 | 8.1147e-03 | -0.13 | 4.7860e-03 | -1.69 | 5.3771e-04 | 0.65 |
| 40 | 5.9738e-03 | 0.87 | 2.2523e-03 | 2.13 | 1.5214e-03 | -2.94 |

Alpha part (the lead field and the direct solve give the same values):

| N | left | order | right | order |
|---|---|---|---|---|
| 10 | 4.4046e-04 |  | 4.4610e-04 |  |
| 14 | 2.2056e-04 | 2.21 | 3.1302e-04 | 1.13 |
| 20 | 4.2193e-05 | 4.71 | 9.9369e-05 | 3.27 |
| 28 | 2.9254e-05 | 1.11 | 9.7180e-05 | 0.07 |
| 40 | 2.0407e-05 | 1.02 | 5.6939e-05 | 1.51 |

- **Lead field and direct solve.** The reciprocity gap is at the level of the
  fixed-point tolerances (at most 2.7e-10 over every group, `layered` and
  `layeredAlpha0`), and the lead field equals the direct solve in every
  row, alpha part included.
- **Total error.** The body-surface ECG converges below second order at the
  intermediate levels (1.82, 1.80) and at second order between `N = 28` and
  `40` (1.95), and the RMS over the surface cells follows it (1.99 at
  `N = 40`). With `-alpha 0` the surface group and both surface field rows
  converge at second order.
- **Local first-order maximum.** With the current, the largest error of the
  field lies in a surface cell at every level (the maxima over all cells and
  over the surface cells coincide), where the layer interface `x = 1.5`
  meets a side face of the box, at `(1.49557 0.2875 0.00541266)` at
  `N = 40`. It converges at about first order, with erratic pairwise orders
  between 0.94 and 1.71. With `-alpha 0` the surface maximum is second
  order, so this local error comes with the current, at the layer interface
  next to the outer boundary. Its cause is not established. A possible
  cause is the explicit non-orthogonal correction across the layer
  conductivity jump at the wall: the faces between the torso layers carry
  the scalar distance-weighted harmonic conductivity, the correction uses
  the standard cell gradient, and the jump-aware gradient of
  `conormalHarmonic` applies only to heart-torso faces.
- **Probes.** The `torso` and `near` probes converge at about second order
  with scatter between levels, consistent with the surface group.
- **Point errors.** The point-evaluated electrode errors converge with
  erratic level-to-level orders in every group, through the `O(h)`
  localisation of each electrode in its cell.
- **Alpha part.** It keeps decreasing to `N = 40`, by a factor of about 20
  over the ladder at `left` and about 8 at `right`, with erratic
  level-to-level orders at both. At `N = 40` it is 4% (`left`) and 12%
  (`right`) of the surface error. No order is measured, and the cause of the
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
  cases are different meshes. The comparison ends at `N = 28`, the last slab
  level. Both ECGs converge to `2 k S`, and their
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
