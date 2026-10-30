# leadFieldConvergence - leadFieldECG

At what order do the lead field and the direct forward solve converge to the
exact potential of the insulated slab at the torso electrodes?

The rows are the static rows of the verifier: the lead field and the direct
solve applied to the exact `Vm` at the heart cell centres, so the error is
that of the passive solve and of the recording alone, not of `Vm_h`. The
electrodes form two groups, all in the torso, where the exact potential is
constant in each slab, so the point and the centroid values agree:

- `torso` (`torsoLeft`, `torsoLeftFar`, `torsoRight`), away from the heart;
- `near` (`nearLeft`, `nearRight`), in torso cells next to the heart.

Each entry is the largest Linf over the samples among the electrodes of the
group. The left slab holds the pin, so the largest error of each group is
that of its right-slab electrode.

Two conductivity configurations, both monodomain with `corrected` snGrad:

- `default` (`-conductivity anisotropic`): `Ge` proportional to `Gi`. The
  heart rows are exact for the exact `Vm`, apart from their coupling across
  the interface, so the error measures the interface, the torso and the
  recording, not the heart operator.
- `unequal` (`-conductivity unequal`): Gi and Ge with different anisotropy
  ratios and the same `xx` values, so the same exact solution. On tetrahedra
  the heart rows carry a truncation error and the passive heart operator is
  exercised as well. On the hexahedra the result equals that of `default`
  (see the case [README](../../../README.md)).

## Ladder

| mesh | N | cells |
|---|---|---|
| hex | 10, 20, 40 | 3000, 24000, 192000 |
| tet (gmsh 4.15.2, seed 1) | 10, 14, 20, 28 | 14792, 39131, 111482, 299957 |

All runs are serial and use `-endTime 0.02` (two samples, `dt = 0.01`). The
lead-field and direct fixed points stop at the `leadFieldFixedPoint`
tolerance of `system/fvSolution`. On hexahedra the lead fields need 5 or 6
passes and each direct solve 2. On tetrahedra the lead fields need 40 to 77
passes and each direct solve 39 to 84.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<configuration>_<mesh>_N<n>`, the names that `../summarize.py` reads. From
this folder:

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
    run default_hex_N$n -mesh hex -N $n
    run unequal_hex_N$n -mesh hex -N $n -conductivity unequal
done
for n in 10 14 20 28
do
    run default_tet_N$n -mesh tet -N $n
    run unequal_tet_N$n -mesh tet -N $n -conductivity unequal
done
python3 ../summarize.py results
```

## Measured behaviour

Errors against the exact potential at the electrode cell centroid, as the
largest Linf over the samples among the electrodes of the group (the exact
potentials reach 1.43; the order is between consecutive levels, with
`h = (cells/3)^(-1/3)`, which is `1/N` on hexahedra).

Hexahedra, `default`:

| N | h | LF torso | order | direct torso | order | LF near | order | direct near | order |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.1 | 1.7585e-02 |  | 1.7585e-02 |  | 1.7585e-02 |  | 1.7585e-02 |  |
| 20 | 0.05 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 |
| 40 | 0.025 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 |

Hexahedra, `unequal`: the same values as `default` to every printed digit.

Tetrahedra, `default`:

| N | h | LF torso | order | direct torso | order | LF near | order | direct near | order |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.05875 | 5.8575e-03 |  | 5.8575e-03 |  | 6.6472e-03 |  | 6.6472e-03 |  |
| 14 | 0.04248 | 3.0903e-03 | 1.97 | 3.0903e-03 | 1.97 | 2.9428e-03 | 2.51 | 2.9428e-03 | 2.51 |
| 20 | 0.02997 | 1.6008e-03 | 1.88 | 1.6008e-03 | 1.88 | 1.6392e-03 | 1.68 | 1.6392e-03 | 1.68 |
| 28 | 0.02155 | 8.1462e-04 | 2.05 | 8.1462e-04 | 2.05 | 7.1188e-04 | 2.53 | 7.1188e-04 | 2.53 |

Tetrahedra, `unequal`:

| N | h | LF torso | order | direct torso | order | LF near | order | direct near | order |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.05875 | 6.9005e-03 |  | 6.9005e-03 |  | 7.5938e-03 |  | 7.5938e-03 |  |
| 14 | 0.04248 | 3.5622e-03 | 2.04 | 3.5622e-03 | 2.04 | 3.4221e-03 | 2.46 | 3.4221e-03 | 2.46 |
| 20 | 0.02997 | 1.8683e-03 | 1.85 | 1.8683e-03 | 1.85 | 1.9919e-03 | 1.55 | 1.9919e-03 | 1.55 |
| 28 | 0.02155 | 9.5185e-04 | 2.04 | 9.5185e-04 | 2.04 | 8.7854e-04 | 2.48 | 8.7854e-04 | 2.48 |

- **Hexahedra.** The lead field and the direct solve agree to every printed
  digit and converge at second order (2.00) in both groups. The ECG between
  the slabs, `torsoRight - torsoLeft`, is second order as well (invariance
  table of [leadFieldOrganJump](../leadFieldOrganJump/README.md)). With
  either configuration the heart rows are exact on this mesh, so these rows
  measure the interface, the torso and the recording.
- **Tetrahedra.** The lead field equals the direct solve to every printed
  digit (gap at most 3.8e-10 relative, see
  [leadFieldReciprocity](../leadFieldReciprocity/README.md)). The centroid
  errors converge at about second order with scatter, 1.85 to 2.05 in the
  torso group: `default` 1.97, 1.88, 2.05, whose heart rows are exact (`Ge`
  proportional to `Gi`); `unequal` 2.04, 1.85, 2.04, with an error about 1.2
  times larger at every level.
- **Electrodes next to the heart.** The `near` group has errors of the size
  of the far torso group and converges at about second order with more
  scatter (1.55 to 2.53).

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
