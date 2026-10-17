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
| hex | 10, 20, 40, 80 | 3000, 24000, 192000, 1536000 |
| tet (gmsh 4.15.2, seed 1) | 10, 14, 20, 28, 40 | 14792, 39131, 111482, 299957, 867138 |

All runs use `-endTime 0.02` (two samples, `dt = 0.01`). On tetrahedra the
lead fields need 36 to 52 passes of their explicit-correction fixed point,
and each direct solve about as many.

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
for n in 10 20 40 80
do
    run default_hex_N$n -mesh hex -N $n
done
for n in 10 20 40
do
    run unequal_hex_N$n -mesh hex -N $n -conductivity unequal
done
for n in 10 14 20 28 40
do
    run default_tet_N$n -mesh tet -N $n
done
for n in 10 14 20 28 40
do
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

| N | LF torso | order | direct torso | order | LF near | order | direct near | order |
|---|---|---|---|---|---|---|---|---|
| 10 | 1.7585e-02 | | 1.7585e-02 | | 1.7585e-02 | | 1.7585e-02 | |
| 20 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 |
| 40 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 |
| 80 | 2.7532e-04 | 2.00 | 2.7532e-04 | 2.00 | 2.7532e-04 | 2.00 | 2.7532e-04 | 2.00 |

Hexahedra, `unequal` (N = 10, 20, 40): the same values as `default` to
1e-12.

Tetrahedra, `default`:

| N | h | LF torso | order | direct torso | order | LF near | order | direct near | order |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.0588 | 6.4987e-03 | | 5.8402e-03 | | 5.0950e-03 | | 6.6177e-03 | |
| 14 | 0.0425 | 3.9020e-03 | 1.57 | 3.0812e-03 | 1.97 | 3.6559e-03 | 1.02 | 2.9405e-03 | 2.50 |
| 20 | 0.0300 | 1.6913e-03 | 2.40 | 1.5963e-03 | 1.88 | 1.8612e-03 | 1.93 | 1.6353e-03 | 1.68 |
| 28 | 0.0216 | 8.0511e-04 | 2.25 | 8.1237e-04 | 2.05 | 6.1093e-04 | 3.38 | 7.1297e-04 | 2.52 |
| 40 | 0.0151 | 3.9453e-04 | 2.02 | 4.0169e-04 | 1.99 | 3.5916e-04 | 1.50 | 4.0338e-04 | 1.61 |

Tetrahedra, `unequal`:

| N | h | LF torso | order | direct torso | order | LF near | order | direct near | order |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.0588 | 4.0509e-03 | | 6.8823e-03 | | 8.5340e-03 | | 7.5662e-03 | |
| 14 | 0.0425 | 6.7166e-04 | 5.54 | 3.5530e-03 | 2.04 | 1.9832e-03 | 4.50 | 3.4220e-03 | 2.45 |
| 20 | 0.0300 | 7.4107e-04 | -0.28 | 1.8637e-03 | 1.85 | 2.0171e-04 | 6.55 | 1.9861e-03 | 1.56 |
| 28 | 0.0216 | 5.3362e-04 | 1.00 | 9.4959e-04 | 2.04 | 8.9572e-04 | -4.52 | 8.7939e-04 | 2.47 |
| 40 | 0.0151 | 3.4775e-04 | 1.21 | 4.7022e-04 | 1.99 | 3.9603e-04 | 2.31 | 4.9056e-04 | 1.65 |

- **Hexahedra.** The lead field and the direct solve agree to every printed
  digit and converge at 2.00 in both groups. The ECG between the slabs,
  `torsoRight - torsoLeft`, has the error of the `LF torso` column. With
  either configuration the heart rows are exact on this mesh, so these
  rows measure the interface, the torso and the recording.
- **Tetrahedra, `default`, direct solve.** It converges at 1.9 to 2.05 in
  the torso group (orders 1.97, 1.88, 2.05, 1.99). Here too the heart rows
  are exact (`Ge` proportional to `Gi`), so the order is that of the
  interface and the torso.
- **Tetrahedra, `unequal`, direct solve.** It converges at 1.85 to 2.04 in
  the torso group (orders 2.04, 1.85, 2.04, 1.99), with an error 10 to 20%
  larger than with `default`, whose heart rows are exact.
- **Tetrahedra, lead field.** With `default` it converges at 1.6 to 2.4 in the
  torso group, 2.02 between `N = 28` and 40, with an error within 30% of the
  direct solve's. With `unequal` its torso error is not monotone (4.1e-3,
  6.7e-4, 7.4e-4, 5.3e-4, 3.5e-4) and no order is measured. The difference
  between the lead field and the direct solve falls from 7.7e-3 to 5.7e-4
  (static rows); see [leadFieldReciprocity](../leadFieldReciprocity/README.md).
- **Electrodes next to the heart.** The `near` group converges as the far
  torso electrodes do, at a less regular order: a torso cell that touches
  the heart does not affect convergence.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
