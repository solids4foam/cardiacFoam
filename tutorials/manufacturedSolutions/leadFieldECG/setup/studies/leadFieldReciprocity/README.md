# leadFieldReciprocity - leadFieldECG

Does the lead-field channel equal the direct forward solve of the same
conductor on the same `Vm`?

The direct solve gives `e^T A^-1 b`, sampled in the electrode cell, with `A`
the forward operator (implicit part and explicit face terms iterated to
their fixed point). The lead field is the discrete adjoint
`z_e = A^-T e`, so the channel `z_e^T b = e^T A^-1 b` equals it for any `A`.
The adjoint matters because `A` is not symmetric on general meshes: the
non-orthogonal correction of the `corrected` snGrad and OpenFOAM's tensor
cross-diffusion `(Sf.G - (Sf.G.n) n) . grad` are explicit terms whose
transpose differs from them, and they vanish only on orthogonal faces whose
normal is an eigenvector of the conductivity (the Cartesian hexahedra with
diagonal conductivities here). `buildLeadFields` transposes them exactly
(probed explicit operator, see its README), and both fixed points stop at the
relative change set by `leadFieldFixedPoint` in `system/fvSolution`
(`tolerance 1e-10`), which sets the size of the gap.

Rows of the verifier summary used here, per electrode (Linf over the
samples):

- `leadField-direct_<electrode>`: the gap on the computed `Vm_h`;
- `leadFieldStatic-directStatic_<electrode>`: the gap on the exact `Vm`.

The summary also has `leadField-direct_relative`, the largest gap over the
electrodes divided by the largest `|direct|`, which the regression checks.

## Ladder

The default configuration is the corrected, anisotropic, monodomain case of
[leadFieldConvergence](../leadFieldConvergence/README.md): hexahedra
`N = 10, 20, 40`, tetrahedra `N = 10, 14, 20, 28`; that study also runs
`unequal` (Gi and Ge with different anisotropy ratios) on the same levels.
This study adds:

| configuration | options | levels |
|---|---|---|
| `isotropic_tet` | `-mesh tet -conductivity isotropic` (non-orthogonal correction only) | 10, 14, 20 |
| `bidomain_hex`, `bidomain_tet` | `-myocardium bidomain` | 10, 20 |
| `parallel_hex`, `parallel_tet` | `parallel` (2 ranks, `numberOfSubdomains 2` in `system/decomposeParDict`) | 20 |

The `parallel` rows run on 2 ranks; every other run is serial.

All runs use `-endTime 0.02` (two samples, `dt = 0.01`). The electrodes are
the case's torso electrodes, in two groups: `torso` (away from the heart)
and `near` (in torso cells next to the heart).

## Execution

Each run is a copy of the case in `results/` (ignored by git). From this
folder:

```bash
mkdir -p results
run()
{
    name=$1; shift
    rsync -a --exclude setup ../../../ results/$name/
    bash results/$name/Allrun "$@" -endTime 0.02
}
for n in 10 14 20
do
    run isotropic_tet_N$n -mesh tet -N $n -conductivity isotropic
done
for n in 10 20
do
    run bidomain_hex_N$n -mesh hex -N $n -myocardium bidomain
    run bidomain_tet_N$n -mesh tet -N $n -myocardium bidomain
done
run parallel_hex_N20 parallel -mesh hex -N 20
run parallel_tet_N20 parallel -mesh tet -N 20
python3 ../summarize.py ../leadFieldConvergence/results results
```

## Measured behaviour

Gap relative to the largest exact potential (1.43), as the largest Linf over
the samples among the electrodes of the group: `dynamic` is
`leadField-direct` on `Vm_h`, `static` is `leadFieldStatic-directStatic` on
the exact `Vm`. Tetrahedra: 14792, 39131, 111482 and 299957 cells at
`N = 10, 14, 20, 28`.

Hexahedra (corrected; the snGrad scheme has no effect there):

| configuration | N | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|
| default | 10 | 2.451e-13 | 2.806e-13 | 7.693e-14 | 3.080e-13 |
|  | 20 | 1.161e-12 | 6.813e-13 | 1.298e-12 | 6.731e-13 |
|  | 40 | 1.481e-12 | 2.955e-12 | 1.248e-12 | 1.988e-12 |
| `unequal` | 10 | 1.637e-13 | 2.584e-13 | 3.400e-13 | 6.903e-13 |
|  | 20 | 6.848e-13 | 6.368e-13 | 1.219e-12 | 1.169e-12 |
|  | 40 | 9.855e-13 | 2.929e-12 | 1.008e-12 | 2.998e-12 |
| bidomain | 10 | 1.470e-13 | 9.779e-14 | 1.206e-13 | 6.607e-14 |
|  | 20 | 2.470e-13 | 8.419e-14 | 2.652e-13 | 9.468e-14 |
| parallel (2 ranks) | 20 | 3.591e-12 | 1.566e-12 | 2.062e-12 | 2.298e-12 |

Tetrahedra:

| configuration | N | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|
| default (corrected, anisotropic) | 10 | 4.646e-11 | 4.142e-11 | 3.756e-11 | 3.455e-11 |
|  | 14 | 1.736e-10 | 1.759e-10 | 2.484e-10 | 2.536e-10 |
|  | 20 | 8.177e-11 | 8.257e-11 | 1.245e-10 | 1.252e-10 |
|  | 28 | 1.129e-10 | 1.137e-10 | 1.102e-10 | 1.099e-10 |
| `unequal` | 10 | 4.822e-11 | 4.305e-11 | 3.212e-11 | 3.344e-11 |
|  | 14 | 3.540e-10 | 3.571e-10 | 2.532e-10 | 2.536e-10 |
|  | 20 | 3.782e-10 | 3.766e-10 | 2.704e-10 | 2.716e-10 |
|  | 28 | 1.006e-10 | 1.008e-10 | 1.105e-10 | 1.105e-10 |
| `isotropic` | 10 | 1.310e-12 | 1.831e-12 | 1.047e-12 | 1.672e-12 |
|  | 14 | 3.337e-12 | 4.007e-12 | 1.783e-12 | 1.935e-12 |
|  | 20 | 2.693e-12 | 2.819e-12 | 3.183e-12 | 3.372e-12 |
| bidomain | 10 | 1.351e-11 | 1.299e-11 | 1.098e-11 | 1.081e-11 |
|  | 20 | 4.039e-11 | 3.970e-11 | 4.418e-11 | 4.488e-11 |
| parallel (2 ranks) | 20 | 8.010e-11 | 8.208e-11 | 1.270e-10 | 1.302e-10 |

- **Hexahedra.** The gap is at round-off, 6.6e-14 to 3.6e-12, and rises
  slightly with the cell count. It is of the same size for `unequal`, for
  the bidomain myocardium and on two ranks.
- **Tetrahedra.** The gap is at the level of the fixed-point tolerances,
  1.1e-11 to 3.8e-10, at every level, for every configuration and both
  electrode groups, with no trend in `N`. With only the non-orthogonal
  correction (`isotropic`) the fixed points take fewer passes (lead fields
  33 to 51, direct solves 28 to 53, against 40 to 77 and 39 to 84 for
  `default` and `unequal`) and the gap is 1.0e-12 to 4.0e-12.
- **Parallel.** The parallel channels agree with the serial ones to 5.8e-12
  (absolute) on both meshes.
- **Electrodes next to the heart.** The `near` group, in torso cells that
  touch the heart, follows the far `torso` group. Electrodes inside the
  heart are not supported (case [README](../../../README.md)).

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
