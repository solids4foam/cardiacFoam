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
`Allrun`); every other run is serial. All runs use the bidomain myocardium
(`-myocardium bidomain`) with `corrected` snGrad, and use
`-endTime 0.02` (two samples, `dt = 0.01`). Each run builds 13 lead fields,
one per electrode. The lead fields and the direct solves use the assembled
matrix through PETSc (GMRES preconditioned by GAMG) and stop on the true
residual of the assembled system at the `residualTolerance 1e-10` of
`system/fvSolution`. On hexahedra the lead fields need 23 to 32 iterations
and each direct solve 24 to 34; on tetrahedra 29 to 36 and 31 to 39. The true
residuals are at most `1.6e-15` (lead fields) and `1.1e-14` (direct
solves).

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
    bash results/$name/Allrun "$@" -myocardium bidomain -endTime 0.02
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

The omnidriver specs for this study are `sweep_hex_layered.json`,
`sweep_hex_slab.json`, `sweep_tet_layered.json` and `sweep_tet_slab.json`, for
the record `manufacturedLeadFieldECG`, which omnidriver does not have. They
name the `Allrun` options as axes (`problem`, `alpha`, `mesh`, `parallel`,
with `myocardium bidomain` in the base;
`-N` is `numberCells` on hexahedra and `tetNumberCells` on tetrahedra) and `-endTime` as the keys
`system/controlDict:endTime` and `writeInterval`.

## Measured behaviour

Errors against the exact potential, as the largest Linf over the samples
among the electrodes of the group, and the error of the direct static field
over the surface cells; the order is between consecutive levels, with
`h = (cells/3)^(-1/3)`, which is `1/N` on hexahedra. Gaps are relative to
the largest exact potential (about 2.05 with `alpha = 0.01`). The lead field
and the direct solve give the same electrode errors to every printed digit,
so the probe and point tables give the lead field only.

### Hexahedra

Reciprocity gap `leadField-direct` (on `Vm_h`) and static, every group, for
`layered` and `layeredAlpha0`: at most 1.2e-12 (N = 10 to 40).

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
| 10 | 5.6510e-14 | 2.6623e-13 | 1.1852e-14 | 3.3751e-14 |
| 20 | 1.0786e-13 | 2.8644e-13 | 1.2768e-13 | 1.8563e-13 |
| 40 | 1.6647e-12 | 1.8150e-12 | 2.8670e-12 | 2.9532e-12 |

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
  (`<= 3.0e-12`) at both probe groups, in the lead field and the direct
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
| 10 | 0.05792 | 3.939e-13 | 3.100e-13 | 6.1177e-03 |  | 6.1177e-03 |  | 3.8674e-03 |  | 1.0075e-02 |  |
| 14 | 0.04238 | 9.916e-13 | 8.005e-13 | 3.0531e-03 | 2.23 | 3.0531e-03 | 2.23 | 2.0339e-03 | 2.06 | 6.2089e-03 | 1.55 |
| 20 | 0.02983 | 1.031e-12 | 9.117e-13 | 1.5876e-03 | 1.86 | 1.5876e-03 | 1.86 | 1.0479e-03 | 1.89 | 4.3815e-03 | 0.99 |
| 28 | 0.02148 | 1.131e-12 | 7.010e-13 | 8.5583e-04 | 1.88 | 8.5583e-04 | 1.88 | 5.6353e-04 | 1.89 | 2.4772e-03 | 1.74 |
| 40 | 0.01508 | 5.425e-12 | 6.159e-12 | 4.1331e-04 | 2.06 | 4.1331e-04 | 2.06 | 2.6832e-04 | 2.10 | 1.7813e-03 | 0.93 |

Body surface with `-alpha 0`:

| N | LF surface | order | `phiStatic_surface` RMS | order | `phiStatic_surface` max | order |
|---|---|---|---|---|---|---|
| 10 | 5.8859e-03 |  | 3.8159e-03 |  | 5.8908e-03 |  |
| 14 | 3.1541e-03 | 2.00 | 2.0512e-03 | 1.99 | 3.2005e-03 | 1.95 |
| 20 | 1.5840e-03 | 1.96 | 1.0354e-03 | 1.95 | 1.5857e-03 | 2.00 |
| 28 | 8.3815e-04 | 1.94 | 5.4643e-04 | 1.94 | 8.4232e-04 | 1.92 |
| 40 | 4.1059e-04 | 2.02 | 2.6798e-04 | 2.01 | 4.1947e-04 | 1.97 |

Probes, `torso` reciprocity gap and lead-field errors:

| N | gap dyn | gap static | LF torso | order | LF near | order |
|---|---|---|---|---|---|---|
| 10 | 1.757e-13 | 1.461e-13 | 6.0582e-03 |  | 6.2693e-03 |  |
| 14 | 8.164e-13 | 6.437e-13 | 3.2297e-03 | 2.01 | 3.1637e-03 | 2.19 |
| 20 | 9.255e-13 | 8.183e-13 | 1.4915e-03 | 2.20 | 1.5657e-03 | 2.00 |
| 28 | 1.174e-12 | 7.487e-13 | 8.0324e-04 | 1.88 | 8.3107e-04 | 1.93 |
| 40 | 5.492e-12 | 4.004e-12 | 4.0038e-04 | 1.97 | 4.0377e-04 | 2.04 |

Static errors at the electrode point (with `-alpha 0` they equal the
centroid errors):

| N | LF surface | order | LF torso | order | LF near | order |
|---|---|---|---|---|---|---|
| 10 | 2.3271e-02 |  | 2.1863e-02 |  | 5.5763e-03 |  |
| 14 | 1.1286e-02 | 2.32 | 1.5397e-02 | 1.12 | 2.0192e-03 | 3.25 |
| 20 | 7.6365e-03 | 1.11 | 2.6488e-03 | 5.01 | 5.6745e-04 | 3.62 |
| 28 | 8.2141e-03 | -0.22 | 4.7857e-03 | -1.80 | 6.3471e-04 | -0.34 |
| 40 | 5.9730e-03 | 0.90 | 2.1857e-03 | 2.22 | 1.5187e-03 | -2.47 |

Alpha part (the lead field and the direct solve give the same values):

| N | left | order | right | order |
|---|---|---|---|---|
| 10 | 4.3298e-04 |  | 1.5523e-04 |  |
| 14 | 2.2088e-04 | 2.15 | 7.7113e-05 | 2.24 |
| 20 | 9.1908e-05 | 2.50 | 9.2567e-05 | -0.52 |
| 28 | 2.9296e-05 | 3.48 | 3.3346e-05 | 3.11 |
| 40 | 2.0728e-05 | 0.98 | 1.0321e-05 | 3.32 |

- **Lead field and direct solve.** The reciprocity gap is at the level of the
  linear-solve tolerances (at most 1.4e-11 relative over every group,
  `layered` and `layeredAlpha0`), and the lead field equals the direct
  solve in every row, alpha part included.
- **Total error.** The body-surface ECG converges below second order at the
  intermediate levels (1.86, 1.88) and at second order between `N = 28` and
  `40` (2.06), and the RMS over the surface cells follows it (2.10 at
  `N = 40`). With `-alpha 0` the surface group and both surface field rows
  converge at second order.
- **Local first-order maximum.** With the current, the largest error of the
  field lies in a surface cell at every level (the maxima over all cells and
  over the surface cells coincide), where the layer interface `x = 1.5`
  meets a side face of the box, at `(1.49557 0.2875 0.00541266)` at
  `N = 40`. It converges at about first order, with erratic pairwise orders
  between 0.93 and 1.74. With `-alpha 0` the surface maximum is second
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
- **Alpha part.** It decreases over the ladder, by a factor of about 20 at
  `left` and about 15 at `right` between `N = 10` and `40`, with erratic
  level-to-level orders at both (at `right` it rises between `N = 14` and
  `20`). At `N = 40` it is 5% (`left`) and 2.5% (`right`) of the surface
  error. No order is measured, and the cause of the
  scatter is not established.

### Invariance: `-alpha 0` against the slab

ECG `torsoRight - torsoLeft`, largest difference over the samples between
`layeredAlpha0` and `slab` at the same `N`, from the channel on `Vm_h` and
from the static channel, and each channel's error against the exact
`2 k S`:

| mesh | N | difference (Vm_h) | difference (static) | error layered | error slab |
|---|---|---|---|---|---|
| hex | 10 | 9.519e-13 | 9.583e-13 | 1.7674e-02 | 1.7674e-02 |
| hex | 20 | 2.369e-13 | 2.389e-13 | 4.5020e-03 | 4.5020e-03 |
| hex | 40 | 4.345e-13 | 4.159e-13 | 1.2024e-03 | 1.2024e-03 |
| tet | 10 | 7.277e-05 | 6.080e-05 | 6.2702e-03 | 6.1974e-03 |
| tet | 14 | 8.359e-05 | 7.300e-05 | 3.4256e-03 | 3.3420e-03 |
| tet | 20 | 1.711e-05 | 1.996e-05 | 1.7859e-03 | 1.8009e-03 |
| tet | 28 | 2.579e-05 | 2.163e-05 | 9.9361e-04 | 9.6782e-04 |

- **Hexahedra.** The two cases have the same cells (five blocks against
  three, the same `h`), and the ECG does not move with the layer
  conductivities: the difference is round-off.
- **Tetrahedra.** The layered mesh has two more internal planes, so the two
  cases are different meshes. The comparison ends at `N = 28`, the last slab
  level. Both ECGs converge to `2 k S`, and their
  difference, 8e-5 to 2e-5, is 1 to 3% of either error.

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
