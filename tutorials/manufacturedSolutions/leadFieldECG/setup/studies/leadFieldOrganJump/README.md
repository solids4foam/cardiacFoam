# leadFieldOrganJump - leadFieldECG

With torso layers of different conductivity and a current across them
(`-problem layered`), at what order do the lead field and the direct forward
solve converge to the exact potential? And is the ECG of the insulated case
(`-alpha 0`) independent of the layer conductivities?

The problem is described in the case [README](../../../README.md): torso
layers `sigma_1 = 0.1` (next to the heart) and `sigma_2 = 0.02` (beyond
`x = -0.5` and `x = 1.5`), a jump of 5, and the surface currents `-alpha` on
`xMin` and `alpha` on `xMax` with `alpha = 0.01`. The torso faces between the
layers take the distance-weighted harmonic mean of the two conductivities, in
the forward operator and therefore in the lead fields.

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
| hex | 10, 20, 40, 80 | 3000, 24000, 192000, 1536000 |
| tet (gmsh 4.15.2, seed 1, `layered_box.geo.template`) | 10, 14, 20, 28, 40 | 15443, 39418, 112979, 302837, 875170 |

All runs are monodomain with `corrected` snGrad and use `-endTime 0.02` (two
samples, `dt = 0.01`).

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
for n in 10 20 40 80
do
    run layered_hex_N$n -problem layered -mesh hex -N $n
    run layeredAlpha0_hex_N$n -problem layered -alpha 0 -mesh hex -N $n
    run slab_hex_N$n -problem slab -mesh hex -N $n
done
for n in 10 14 20 28 40
do
    run layered_tet_N$n -problem layered -mesh tet -N $n
    run layeredAlpha0_tet_N$n -problem layered -alpha 0 -mesh tet -N $n
    run slab_tet_N$n -problem slab -mesh tet -N $n
done
# Control without explicit face terms (see the tetrahedra below)
for n in 10 14 20
do
    run orthIsoLayered_tet_N$n -problem layered -mesh tet -N $n \
        -snGrad orthogonal -conductivity isotropic
    run orthIsoLayeredAlpha0_tet_N$n -problem layered -alpha 0 -mesh tet -N $n \
        -snGrad orthogonal -conductivity isotropic
done
python3 ../summarize.py results
```

## Measured behaviour

Errors against the exact potential, as the largest Linf over the samples
among the electrodes of the group; the order is between consecutive levels,
with `h = (cells/3)^(-1/3)`, which is `1/N` on hexahedra. Gaps are relative to
the largest exact potential (1.79 with `alpha = 0.01`).

### Hexahedra

Reciprocity gap `leadField-direct` (on `Vm_h`) and static, both groups: at
most `6.7e-12` (N = 10 to 80); the lead field and the direct solve agree to
every printed digit below.

Static errors at the cell centroid (`torso` and `near` give the same values,
set by the right-slab electrodes):

| N | LF torso | order | direct torso | order | `-alpha 0` torso | order |
|---|---|---|---|---|---|---|
| 10 | 1.3098e-02 | | 1.3098e-02 | | 1.7585e-02 | |
| 20 | 2.1598e-03 | 2.60 | 2.1598e-03 | 2.60 | 4.4029e-03 | 2.00 |
| 40 | 2.5801e-05 | 6.39 | 2.5801e-05 | 6.39 | 1.1012e-03 | 2.00 |
| 80 | 2.8681e-04 | -3.47 | 2.8681e-04 | -3.47 | 2.7532e-04 | 2.00 |

Alpha part (lead field; the direct solve gives the same values):

| N | left | right | order (right) |
|---|---|---|---|
| 10 | 4.3e-14 | 4.4862e-03 | |
| 20 | 1.2e-13 | 2.2431e-03 | 1.00 |
| 40 | 1.7e-12 | 1.1215e-03 | 1.00 |
| 80 | 4.3e-11 | 5.6077e-04 | 1.00 |

Static errors at the electrode point (`torso`): 1.0650e-02, 1.1310e-02,
4.4000e-03, 1.2750e-03. With `alpha != 0` the potential has a slope of up to
0.5 in the outer layers, so the point error is the localisation
`phi(x_e) - phi(x_c)` of the electrode in its cell, `O(h)` and not monotone.

- **The layer interfaces are exact on the hexahedra.** The `left` alpha part,
  whose current crosses `x = -0.5` from the pinned layer, is round-off
  (`<= 4.3e-11`) in the lead field and the direct solve: the harmonic face
  conductivity gives the exact flux of a piecewise linear potential across a
  jump of 5.
- **The heart-torso interface is first order on the hexahedra.** The `right`
  alpha part converges at 1.00. It equals, to ten digits,
  `alpha h (1/Ge_xx - 1/(Gi + Ge)_xx)` (`4.4861838189e-03` at `N = 10`
  against `4.4861838189e-03` predicted): the potential offset of the two
  heart-torso faces, which carry the extracellular conductivity `Ge` only
  (the intracellular part of the operator is insulated there), while the
  exact potential has the slope `alpha/(Gi + Ge)_xx` in the heart and so the
  intracellular flux `alpha Gi_xx/(Gi + Ge)_xx` at `x = 0` and `x = 1`.
- **The total error has no clean order.** The `O(h^2)` slab part (the
  `-alpha 0` column) and the `O(h)` alpha part have opposite signs at the
  right electrodes; they nearly cancel at `N = 40`, and the alpha part
  dominates from `N = 80`.

### Tetrahedra

| N | h | gap dyn | gap static | LF torso | order | direct torso | order | LF near | direct near |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.0579 | 2.3e-03 | 1.6e-03 | 6.6405e-03 | | 3.7456e-03 | | 3.9595e-03 | 4.0063e-03 |
| 14 | 0.0424 | 3.8e-04 | 6.7e-04 | 3.7113e-04 | 9.23 | 1.5545e-03 | 2.82 | 8.6397e-04 | 1.4928e-03 |
| 20 | 0.0298 | 8.9e-04 | 8.8e-04 | 1.5633e-03 | -4.10 | 2.8059e-04 | 4.88 | 1.6156e-03 | 3.2376e-04 |
| 28 | 0.0215 | 1.3e-03 | 9.3e-04 | 1.6576e-03 | -0.18 | 7.8267e-05 | 3.88 | 1.1822e-03 | 4.0400e-05 |
| 40 | 0.0151 | 3.0e-04 | 7.4e-05 | 3.2245e-04 | 4.63 | 2.2137e-04 | -2.94 | 3.2099e-04 | 2.4892e-04 |

(`gap`: the `torso` reciprocity gap on `Vm_h` and on the exact `Vm`.)

Static errors at the electrode point (`torso`), lead field: 2.1278e-02,
1.5463e-02, 2.1713e-03, 3.1409e-03, 1.4645e-03; direct: 2.4175e-02,
1.5393e-02, 2.6485e-03, 4.7848e-03, 1.5659e-03.

Alpha part:

| N | LF left | LF right | direct left | order | direct right | order |
|---|---|---|---|---|---|---|
| 10 | 2.1664e-03 | 1.2306e-03 | 4.6069e-04 | | 2.4310e-03 | |
| 14 | 8.2195e-04 | 2.4316e-03 | 2.2832e-04 | 2.25 | 1.6823e-03 | 1.18 |
| 20 | 1.5679e-03 | 4.6566e-04 | 4.7255e-05 | 4.49 | 1.2861e-03 | 0.77 |
| 28 | 1.1573e-03 | 8.6693e-04 | 3.9780e-05 | 0.52 | 9.0166e-04 | 1.08 |
| 40 | 7.3169e-05 | 6.8372e-04 | 1.9220e-05 | 2.06 | 6.5066e-04 | 0.92 |

`-snGrad orthogonal -conductivity isotropic` (no explicit face terms,
`N = 10, 14, 20`): the reciprocity gap is at most `1.2e-12`, and the lead
field equals the direct solve in every row, alpha part included.

- **Direct solve, layer interface.** The `left` alpha part falls from 4.6e-4
  to 1.9e-5 (orders 2.25, 4.49, 0.52, 2.06): it converges, at an irregular
  order of about 2.
- **Direct solve, heart-torso interface.** The `right` alpha part converges at
  about 1 (1.18, 0.77, 1.08, 0.92), as on the hexahedra.
- **Lead field.** Its alpha part does not converge from `N = 10` to 28
  (left 2.2e-3, 8.2e-4, 1.6e-3, 1.2e-3) and drops to 7.3e-5 at `N = 40`;
  no order is measured. The reciprocity gap behaves the same way (about 1e-3
  up to `N = 28`, 3.0e-4 at 40). Without the explicit face terms the gap is
  round-off, so it is the asymmetry of the non-orthogonal and cross-diffusion
  corrections: the lead field is the transpose of the corrected forward
  operator, and here the torso potential has a gradient on which those
  corrections act. In the insulated slab the torso potential is constant and
  this does not show (`-alpha 0`: gap 5.6e-4 down to 3.9e-4, static 3.0e-4
  down to 5.3e-5).
- The total static error is not monotone for the lead field, and for the
  direct solve it mixes the `O(h)` heart-interface part with the rest; no
  order is measured for either.

### Invariance: `-alpha 0` against the slab

ECG `torsoRight - torsoLeft`, largest difference over the samples between
`layeredAlpha0` and `slab` at the same `N`, from the channel on `Vm_h` and
from the static channel, and each channel's error against the exact
`2 k S`:

| mesh | N | difference (Vm_h) | difference (static) | error layered | error slab |
|---|---|---|---|---|---|
| hex | 10 | 3.7e-13 | 3.5e-13 | 1.7645e-02 | 1.7645e-02 |
| hex | 20 | 1.3e-12 | 1.3e-12 | 4.4945e-03 | 4.4945e-03 |
| hex | 40 | 1.7e-12 | 1.7e-12 | 1.2006e-03 | 1.2006e-03 |
| hex | 80 | 4.5e-12 | 4.5e-12 | 3.7668e-04 | 3.7668e-04 |
| tet | 10 | 2.4e-03 | 1.0e-03 | 7.3314e-03 | 9.7193e-03 |
| tet | 14 | 2.1e-03 | 1.1e-03 | 4.1766e-03 | 6.2965e-03 |
| tet | 20 | 1.1e-03 | 4.2e-04 | 2.4200e-03 | 3.4851e-03 |
| tet | 28 | 4.7e-04 | 1.5e-04 | 1.5889e-03 | 2.0574e-03 |
| tet | 40 | 2.5e-04 | 7.4e-05 | 1.1011e-03 | 1.3496e-03 |

- **Hexahedra.** The two cases have the same cells (five blocks against
  three, the same `h`), and the ECG does not move with the layer
  conductivities: the difference is round-off.
- **Tetrahedra.** The layered mesh has two more internal planes, so the two
  cases are different meshes. Both ECGs converge to `2 k S`, and their
  difference falls with them (2.4e-3 to 2.5e-4).

### Linear interpolation of the layer conductivity

The torso face rule has no switch: torso faces between different
conductivities always take the distance-weighted harmonic mean, so a linear
control was not run. On the hexahedra the harmonic mean is exact for the
piecewise linear potential (the `left` rows). With the linear mean
`(sigma_1 + sigma_2)/2` instead, the face between the layers would carry the
flux of a resistance `h/(sigma_1 + sigma_2) * 2` instead of
`h/2 (1/sigma_1 + 1/sigma_2)`, an offset of
`alpha h ((1/sigma_1 + 1/sigma_2)/2 - 2/(sigma_1 + sigma_2))` per interface
(`13.3 alpha h` here): first order at the jump.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
