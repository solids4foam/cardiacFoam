# leadFieldReciprocity - leadFieldECG

Does the lead-field channel equal the direct forward solve of the same
conductor on the same `Vm`?

The lead field `z_e` solves `A z_e = delta_e` with the converged forward
operator `A` (implicit part and explicit face terms iterated to their fixed
point), and the channel is `z_e^T b = e^T A^-T b`. The direct solve gives
`e^T A^-1 b`, sampled in the electrode cell. The two are equal when `A` is
symmetric. That holds when the explicit face terms of the operator vanish:

- the non-orthogonal correction of the `corrected` snGrad, zero on
  orthogonal faces;
- OpenFOAM's tensor cross-diffusion `(Sf.G - (Sf.G.n) n) . grad`, which the
  tensor Laplacian adds whatever the snGrad scheme. It is zero when the face
  normal is an eigenvector of the conductivity.

On the Cartesian hexahedra with diagonal conductivities both terms vanish.
On tetrahedra the non-orthogonal correction is removed by
`-snGrad orthogonal`, and the cross-diffusion only by an isotropic heart
conductivity (`-conductivity isotropic`). Otherwise the
channel is the solution of `A^T`, whose explicit part is not a consistent
discretisation in general.

`-snGrad orthogonal` is itself an inconsistent scheme on tetrahedra: it
takes the face-normal gradient as the difference of the two cell values
over their distance, which is not the normal gradient where the line between
the centres is not along the face normal. Only `orthogonal` together with
`isotropic` leaves no explicit term, so that configuration is the symmetry
control of the machinery (lead field, recording, transpose). The
`orthogonal` rows with the anisotropic heart combine an inconsistent
symmetric operator with the transposed cross-diffusion and do not isolate
either.

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
| `orthogonal_tet` | `-mesh tet -snGrad orthogonal` | 10, 14, 20, 28 |
| `isotropic_tet` | `-mesh tet -conductivity isotropic` | 10, 14, 20 |
| `orthogonalIsotropic_tet` | `-mesh tet -snGrad orthogonal -conductivity isotropic` | 10, 14, 20, 28 |
| `bidomain_hex`, `bidomain_tet` | `-myocardium bidomain` | 10, 20 |
| `parallel_hex` | `parallel` (2 ranks) | 20 |

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
for n in 10 14 20 28
do
    run orthogonal_tet_N$n -mesh tet -N $n -snGrad orthogonal
    run orthogonalIsotropic_tet_N$n -mesh tet -N $n -snGrad orthogonal -conductivity isotropic
done
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
| default | 10 | 7.524e-14 | 1.467e-13 | 2.101e-13 | 3.183e-13 |
| | 20 | 5.273e-13 | 1.573e-12 | 1.457e-13 | 1.192e-12 |
| | 40 | 6.606e-13 | 1.424e-12 | 1.549e-12 | 1.832e-12 |
| | 80 | 2.403e-12 | 2.402e-12 | 2.867e-12 | 1.857e-12 |
| `unequal` | 10 | 1.767e-13 | 2.621e-13 | 3.129e-13 | 2.254e-13 |
| | 20 | 1.820e-12 | 1.268e-12 | 7.050e-13 | 5.582e-13 |
| | 40 | 1.572e-12 | 2.853e-12 | 1.386e-12 | 1.799e-12 |
| bidomain | 10 | 4.804e-14 | 6.172e-14 | 5.087e-14 | 2.624e-14 |
| | 20 | 1.883e-13 | 3.431e-13 | 1.051e-13 | 2.105e-13 |
| parallel (2 ranks) | 20 | 3.250e-12 | 2.585e-12 | 2.575e-12 | 3.516e-12 |

Tetrahedra:

| configuration | N | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|
| default (corrected, anisotropic) | 10 | 2.243e-03 | 4.611e-04 | 3.209e-03 | 3.276e-03 |
| | 14 | 1.905e-03 | 5.747e-04 | 1.748e-03 | 5.008e-04 |
| | 20 | 1.085e-03 | 6.653e-05 | 1.080e-03 | 1.582e-04 |
| | 28 | 7.269e-04 | 5.083e-06 | 7.160e-04 | 7.144e-05 |
| | 40 | 5.678e-04 | 5.015e-06 | 5.141e-04 | 3.096e-05 |
| `unequal` (corrected) | 10 | 5.543e-03 | 7.655e-03 | 8.785e-03 | 1.100e-02 |
| | 14 | 1.436e-03 | 2.958e-03 | 2.084e-03 | 3.454e-03 |
| | 20 | 6.034e-04 | 1.824e-03 | 4.812e-04 | 1.532e-03 |
| | 28 | 1.936e-04 | 1.038e-03 | 3.989e-04 | 1.243e-03 |
| | 40 | 3.198e-04 | 5.727e-04 | 2.421e-04 | 6.207e-04 |
| `orthogonal` (anisotropic) | 10 | 1.321e-03 | 4.300e-03 | 3.473e-03 | 3.853e-03 |
| | 14 | 4.301e-03 | 9.662e-03 | 4.378e-03 | 9.775e-03 |
| | 20 | 3.550e-03 | 9.944e-03 | 4.085e-03 | 1.064e-02 |
| | 28 | 3.052e-03 | 1.061e-02 | 3.321e-03 | 1.094e-02 |
| `isotropic` (corrected) | 10 | 6.530e-03 | 7.677e-03 | 1.149e-02 | 1.271e-02 |
| | 14 | 2.830e-03 | 3.677e-03 | 2.793e-03 | 3.578e-03 |
| | 20 | 1.555e-03 | 2.190e-03 | 1.892e-03 | 2.496e-03 |
| `orthogonalIsotropic` | 10 | 6.820e-13 | 9.885e-13 | 6.048e-13 | 7.472e-13 |
| | 14 | 4.339e-13 | 5.259e-13 | 3.865e-13 | 7.616e-13 |
| | 20 | 7.407e-13 | 6.959e-13 | 7.851e-13 | 7.762e-13 |
| | 28 | 2.726e-12 | 2.374e-12 | 1.392e-12 | 1.362e-12 |
| bidomain (corrected, anisotropic) | 10 | 1.601e-04 | 3.729e-04 | 1.482e-03 | 1.498e-03 |
| | 20 | 3.383e-04 | 1.916e-04 | 2.704e-04 | 2.330e-04 |

Lead field against the exact potential at the centroid (largest Linf in the
group) for the tetrahedral controls, with the direct solve for comparison:

| configuration | N | LF torso | direct torso | LF near | direct near |
|---|---|---|---|---|---|
| `orthogonal` | 10 | 9.8287e-04 | 5.1591e-03 | 4.6258e-03 | 5.6628e-03 |
| | 14 | 1.1115e-02 | 2.6856e-03 | 1.1368e-02 | 2.5943e-03 |
| | 20 | 1.2832e-02 | 1.3711e-03 | 1.3806e-02 | 1.3868e-03 |
| | 28 | 1.4453e-02 | 7.0686e-04 | 1.4997e-02 | 6.3427e-04 |
| `isotropic` | 10 | 5.3137e-03 | 5.6511e-03 | 1.1696e-02 | 6.4581e-03 |
| | 14 | 2.2724e-03 | 2.9801e-03 | 2.2614e-03 | 2.8491e-03 |
| | 20 | 1.5832e-03 | 1.5451e-03 | 1.9734e-03 | 1.5914e-03 |
| `orthogonalIsotropic` | 10 | 5.0155e-03 | 5.0155e-03 | 5.5172e-03 | 5.5172e-03 |
| | 14 | 2.6073e-03 | 2.6073e-03 | 2.5283e-03 | 2.5283e-03 |
| | 20 | 1.3328e-03 | 1.3328e-03 | 1.3504e-03 | 1.3504e-03 |
| | 28 | 6.8719e-04 | 6.8719e-04 | 6.1319e-04 | 6.1319e-04 |

- **Hexahedra.** The gap is at the linear-solver level, 2.6e-14 to 3.5e-12,
  and grows slowly with the cell count, as the round-off of the 1e-15
  absolute PCG tolerance does. It is the same for `unequal`, for the
  bidomain myocardium and on two ranks. The parallel channels agree with the
  serial ones to 4.6e-12 (absolute).
- **Tetrahedra, `orthogonalIsotropic`.** With no explicit face term the gap
  is at the solver level (below 2.8e-12), and the lead field equals the
  direct solve to every printed digit. This is the symmetry control of the
  machinery. The order of about 2 in this problem comes from the slab
  (heart rows exact for `Ge` proportional to `Gi`, exact torso potential
  constant); it does not make `orthogonal` a consistent scheme on
  tetrahedra.
- **Tetrahedra, `orthogonal` with anisotropic conductivity.** The scheme is
  inconsistent on these meshes and the tensor cross-diffusion of the
  anisotropic heart stays explicit, so `A` is not symmetric. The gap stays
  between 1.3e-3 and 1.1e-2 and the lead field moves away from the exact
  potential (9.8e-4 to 1.5e-2). The direct solve's order 2 here follows from
  the slab, as above. These rows do not separate the inconsistent symmetric
  operator from the transposed cross-diffusion.
- **Tetrahedra, `isotropic` with the non-orthogonal correction.** The gap
  falls (6.5e-3 to 1.6e-3 dynamic, 7.7e-3 to 2.2e-3 static in the torso
  group), and the lead field converges with the direct solve.
- **Tetrahedra, default and `unequal` (both explicit terms).** The gap falls
  in both groups, irregularly: default torso dynamic 2.2e-3 to 5.7e-4 and
  static 4.6e-4 to 5.0e-6 from `N = 10` to 40; `unequal`
  torso dynamic 5.5e-3 to 3.2e-4 and static 7.7e-3 to 5.7e-4 from `N = 10`
  to 40.
- **Tetrahedra, bidomain.** Two levels only; the gap is between 1.6e-4 and
  1.5e-3, the same size as the monodomain default's.
- **Electrodes next to the heart.** The `near` group, in torso cells that
  touch the heart, follows the far `torso` group in every configuration.
  Electrodes inside the heart are not supported (case
  [README](../../../README.md)).

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
