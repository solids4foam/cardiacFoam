# leadFieldConvergence - leadFieldECG

At what order do the lead field and the direct forward solve converge to the
exact potential of the insulated slab on the body surface?

The rows are the static rows of the verifier: the lead field and the direct
solve applied to the exact `Vm` at the heart cell centres, so the error is
that of the passive solve and of the recording alone, not of `Vm_h`. Every
electrode lies in the torso, where the exact potential is constant in each
slab, so the point and the centroid values agree. The body-surface
electrodes are the ECG and carry the convergence result; the interior
electrodes are probes that show where an error comes from:

- `surface` (`surface1` to `surface8`), the body-surface electrodes, 0.003
  inside the outer faces of the torso slabs;
- `torso` (`torsoLeft`, `torsoLeftFar`, `torsoRight`), interior probes away
  from the heart;
- `near` (`nearLeft`, `nearRight`), interior probes in torso cells next to
  the heart.

Each electrode entry is the largest Linf over the samples among the
electrodes of the group. A single electrode can sit near a zero of the
error and give erratic orders between levels; the group maximum does not
depend on one point. The field row `phiStatic_surface` is the RMS and the
maximum error of the direct static potential over the cells of the outer
boundary. The left slab holds the pin, so the largest error of each
electrode group is that of a right-slab electrode.

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

All runs are serial and use `-endTime 0.02` (two samples, `dt = 0.01`).
Each run builds 13 lead fields, one per electrode. The lead-field and direct
fixed points stop at the `leadFieldFixedPoint` tolerance of
`system/fvSolution`. On hexahedra the lead fields need 5 or 6 passes and
each direct solve 2. On tetrahedra the lead fields need 40 to 77 passes and
each direct solve 39 to 84.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<configuration>_<mesh>_N<n>`, the names that `../summarize.py` reads. From
this folder:

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
potentials reach 1.43), and the error of the direct static field over the
surface cells; the order is between consecutive levels, with
`h = (cells/3)^(-1/3)`, which is `1/N` on hexahedra. The lead field and the
direct solve agree to every printed digit in every row, so the probe tables
give the lead field only.

### Body surface

Hexahedra, `default`:

| N | h | LF surface | order | direct surface | order | `phiStatic_surface` RMS | order | `phiStatic_surface` max | order |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.1 | 1.7585e-02 |  | 1.7585e-02 |  | 1.1470e-02 |  | 1.7585e-02 |  |
| 20 | 0.05 | 4.4029e-03 | 2.00 | 4.4029e-03 | 2.00 | 2.8772e-03 | 2.00 | 4.4029e-03 | 2.00 |
| 40 | 0.025 | 1.1012e-03 | 2.00 | 1.1012e-03 | 2.00 | 7.2023e-04 | 2.00 | 1.1012e-03 | 2.00 |

Hexahedra, `unequal`: the same values as `default` to every printed digit.

Tetrahedra, `default`:

| N | h | LF surface | order | direct surface | order | `phiStatic_surface` RMS | order | `phiStatic_surface` max | order |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.05875 | 5.8315e-03 |  | 5.8315e-03 |  | 3.7627e-03 |  | 5.8465e-03 |  |
| 14 | 0.04248 | 3.0830e-03 | 1.97 | 3.0830e-03 | 1.97 | 2.0003e-03 | 1.95 | 3.1204e-03 | 1.94 |
| 20 | 0.02997 | 1.5989e-03 | 1.88 | 1.5989e-03 | 1.88 | 1.0391e-03 | 1.88 | 1.6101e-03 | 1.90 |
| 28 | 0.02155 | 8.1598e-04 | 2.04 | 8.1598e-04 | 2.04 | 5.2882e-04 | 2.05 | 8.2201e-04 | 2.04 |

Tetrahedra, `unequal`:

| N | h | LF surface | order | direct surface | order | `phiStatic_surface` RMS | order | `phiStatic_surface` max | order |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.05875 | 6.8607e-03 |  | 6.8607e-03 |  | 4.4114e-03 |  | 6.9362e-03 |  |
| 14 | 0.04248 | 3.5545e-03 | 2.03 | 3.5545e-03 | 2.03 | 2.3011e-03 | 2.01 | 3.5742e-03 | 2.04 |
| 20 | 0.02997 | 1.8689e-03 | 1.84 | 1.8689e-03 | 1.84 | 1.2104e-03 | 1.84 | 1.8815e-03 | 1.84 |
| 28 | 0.02155 | 9.5308e-04 | 2.04 | 9.5308e-04 | 2.04 | 6.1726e-04 | 2.04 | 9.6692e-04 | 2.02 |

### Interior probes

Lead field, hexahedra (`default` and `unequal`): `torso` and `near` have the
values and orders of the surface group to every printed digit.

Lead field, tetrahedra:

| configuration | N | LF torso | order | LF near | order |
|---|---|---|---|---|---|
| `default` | 10 | 5.8575e-03 |  | 6.6472e-03 |  |
|  | 14 | 3.0903e-03 | 1.97 | 2.9428e-03 | 2.51 |
|  | 20 | 1.6008e-03 | 1.88 | 1.6392e-03 | 1.68 |
|  | 28 | 8.1462e-04 | 2.05 | 7.1188e-04 | 2.53 |
| `unequal` | 10 | 6.9005e-03 |  | 7.5938e-03 |  |
|  | 14 | 3.5622e-03 | 2.04 | 3.4221e-03 | 2.46 |
|  | 20 | 1.8683e-03 | 1.85 | 1.9919e-03 | 1.55 |
|  | 28 | 9.5185e-04 | 2.04 | 8.7854e-04 | 2.48 |

### Reading

- **Hexahedra.** The body-surface ECG converges at second order in the
  lead field and the direct solve, and the error over the surface cells
  converges at second order in the RMS and the maximum. The lead field and
  the direct solve agree to every printed digit. Every torso cell of a slab
  carries the same offset, so the probes give the surface values. With
  either configuration the heart rows are exact on this mesh, so these rows
  measure the interface, the torso and the recording. The ECG between the
  slabs, `torsoRight - torsoLeft`, is second order as well (invariance
  table of [leadFieldOrganJump](../leadFieldOrganJump/README.md)).
- **Tetrahedra.** The body-surface ECG converges at second order with
  scatter between levels, in both configurations, and the RMS and maximum
  over the surface cells follow the same level-to-level pattern, so the
  scatter belongs to the mesh ladder and not to the electrode positions.
  `unequal`, whose heart rows carry a truncation error, has an error about
  1.2 times that of `default` at every level, with the same orders. The lead
  field equals the direct solve to every printed digit (gap at most
  5.9e-10 relative, see [leadFieldReciprocity](../leadFieldReciprocity/README.md)).
- **Probes.** The `torso` probes have errors of the size of the surface
  group and the same orders. The `near` probes, next to the heart, have
  errors of the same size and converge at about second order with larger
  scatter between levels: a single electrode next to the interface does not
  give a steady pairwise order, while the surface group and the surface
  field row do.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
