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
(probed explicit operator, see its README), and both fixed points stop at a
relative change of `1e-10`, which sets the size of the gap.

Rows of the verifier summary used here, per electrode (Linf over the
samples):

- `leadField-direct_<electrode>`: the gap on the computed `Vm_h`;
- `leadFieldStatic-directStatic_<electrode>`: the gap on the exact `Vm`.

The summary also has `leadField-direct_relative`, the largest gap over the
electrodes divided by the largest `|direct|`, which the regression checks.

## Ladder

The default configuration is the corrected, anisotropic, monodomain case of
[leadFieldConvergence](../leadFieldConvergence/README.md): hexahedra
`N = 10, 20, 40, 80`, tetrahedra `N = 10, 14, 20, 28, 40`; that study also
runs `unequal` (Gi and Ge with different anisotropy ratios) on hexahedra
`N = 10, 20, 40` and tetrahedra `N = 10, 14, 20, 28, 40`. This study adds:

| configuration | options | levels |
|---|---|---|
| `isotropic_tet` | `-mesh tet -conductivity isotropic` (non-orthogonal correction only) | 10, 14, 20 |
| `bidomain_hex`, `bidomain_tet` | `-myocardium bidomain` | 10, 20 |
| `parallel_hex`, `parallel_tet` | `parallel` (2 ranks) | 20 |

The `N = 40` tetrahedra and the `N = 40` and `80` hexahedra were run on 8
ranks (`parallel` with `numberOfSubdomains 8`); the other levels serially,
apart from the `parallel` rows.

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
the exact `Vm`. Tetrahedra: 14792, 39131, 111482, 299957 and 867138 cells at
`N = 10, 14, 20, 28, 40`.

Hexahedra (corrected; the snGrad scheme has no effect there):

| configuration | N | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|
| default | 10 | 5.201e-13 | 3.563e-13 | 4.127e-13 | 5.534e-14 |
| | 20 | 6.041e-13 | 7.734e-13 | 1.147e-13 | 1.820e-13 |
| | 40 | 4.860e-12 | 5.588e-12 | 1.383e-12 | 3.459e-12 |
| | 80 | 1.812e-11 | 2.451e-11 | 2.420e-11 | 1.580e-11 |
| `unequal` | 10 | 2.001e-13 | 1.584e-13 | 2.504e-13 | 2.038e-13 |
| | 20 | 1.248e-12 | 4.972e-13 | 9.975e-13 | 6.416e-13 |
| | 40 | 5.135e-12 | 6.192e-12 | 4.603e-12 | 4.010e-12 |
| bidomain | 10 | 3.321e-14 | 4.828e-14 | 2.301e-14 | 3.731e-14 |
| | 20 | 3.382e-13 | 2.772e-13 | 2.705e-13 | 4.354e-13 |
| parallel (2 ranks) | 20 | 2.926e-12 | 2.127e-12 | 1.788e-12 | 9.756e-13 |

Tetrahedra:

| configuration | N | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|
| default (corrected, anisotropic) | 10 | 3.758e-11 | 3.355e-11 | 3.307e-11 | 3.061e-11 |
| | 14 | 9.158e-11 | 9.317e-11 | 2.001e-10 | 2.005e-10 |
| | 20 | 2.138e-11 | 1.842e-11 | 4.580e-11 | 4.599e-11 |
| | 28 | 9.496e-11 | 9.580e-11 | 1.126e-10 | 1.132e-10 |
| | 40 | 7.475e-12 | 6.129e-12 | 1.258e-11 | 1.086e-11 |
| `unequal` | 10 | 3.858e-11 | 3.399e-11 | 3.371e-11 | 3.201e-11 |
| | 14 | 1.400e-10 | 1.345e-10 | 2.368e-10 | 2.286e-10 |
| | 20 | 7.681e-11 | 7.763e-11 | 2.832e-10 | 2.842e-10 |
| | 28 | 1.935e-10 | 1.934e-10 | 4.887e-11 | 5.051e-11 |
| | 40 | 7.779e-11 | 7.667e-11 | 5.093e-11 | 5.252e-11 |
| `isotropic` | 10 | 1.444e-12 | 2.553e-12 | 1.023e-12 | 2.157e-12 |
| | 14 | 1.498e-12 | 2.124e-12 | 1.354e-12 | 6.959e-13 |
| | 20 | 2.401e-12 | 3.304e-12 | 2.292e-12 | 2.498e-12 |
| bidomain | 10 | 1.140e-11 | 1.151e-11 | 9.249e-12 | 8.861e-12 |
| | 20 | 2.107e-11 | 2.049e-11 | 3.306e-11 | 3.381e-11 |
| parallel (2 ranks) | 20 | 2.016e-11 | 2.181e-11 | 4.565e-11 | 4.837e-11 |

- **Hexahedra.** The gap is at the linear-solver level, 2.3e-14 to 2.5e-11,
  and grows slowly with the cell count, as the round-off of the 1e-15
  absolute PCG tolerance does. It is the same for `unequal`, for the
  bidomain myocardium and on two ranks.
- **Tetrahedra.** The gap is at the level of the two `1e-10` fixed-point
  tolerances, 6e-12 to 2.8e-10, at every level, for every configuration and
  both electrode groups, with no trend in `N`. With only the non-orthogonal
  correction (`isotropic`) the fixed points converge faster and the gap is
  about 2e-12.
- **Parallel.** The parallel channels agree with the serial ones to 3.8e-12
  (absolute) on both meshes.
- **Electrodes next to the heart.** The `near` group, in torso cells that
  touch the heart, follows the far `torso` group. Electrodes inside the
  heart are not supported (case [README](../../../README.md)).

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
