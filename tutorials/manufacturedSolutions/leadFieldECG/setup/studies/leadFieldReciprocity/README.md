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
| `bidomain_hex` | `-mesh hex -myocardium bidomain` | 10, 20, 40 |
| `bidomain_tet` | `-mesh tet -myocardium bidomain` | 10, 14, 20, 28 |
| `parallel_hex`, `parallel_tet` | `parallel` (2 ranks, `numberOfSubdomains 2` in `system/decomposeParDict`) | 20 |
| `fixedPointTight_tet` | `-mesh tet -fixedPointTolerance 1e-12` | 20 |

The `parallel` rows run on 2 ranks; every other run is serial.

All runs use `-endTime 0.02` (two samples, `dt = 0.01`) and build 13 lead
fields, one per electrode. The electrodes are those of the case, in three
groups: `surface` (`surface1` to `surface8`, the body-surface electrodes of
the ECG), and the interior probes `torso` (away from the heart) and `near`
(in torso cells next to the heart).

## Execution

Each run is a copy of the case in `results/` (ignored by git). From this
folder:

```bash
mkdir -p results
run()
{
    name=$1; shift
    rsync -a --exclude setup --exclude regression --exclude sphericalShell \
        ../../../ results/$name/
    bash results/$name/Allrun "$@" -endTime 0.02
}
for n in 10 14 20
do
    run isotropic_tet_N$n -mesh tet -N $n -conductivity isotropic
done
for n in 10 20 40
do
    run bidomain_hex_N$n -mesh hex -N $n -myocardium bidomain
done
for n in 10 14 20 28
do
    run bidomain_tet_N$n -mesh tet -N $n -myocardium bidomain
done
run parallel_hex_N20 parallel -mesh hex -N 20
run parallel_tet_N20 parallel -mesh tet -N 20
run fixedPointTight_tet_N20 -mesh tet -N 20 -fixedPointTolerance 1e-12
python3 ../summarize.py ../leadFieldConvergence/results results
```

## Measured behaviour

Gap relative to the largest exact potential (1.43), as the largest Linf over
the samples among the electrodes of the group: `dynamic` is
`leadField-direct` on `Vm_h`, `static` is `leadFieldStatic-directStatic` on
the exact `Vm`. Tetrahedra: 14792, 39131, 111482 and 299957 cells at
`N = 10, 14, 20, 28`.

Hexahedra (corrected; the snGrad scheme has no effect there):

| configuration | N | surface dynamic | surface static | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|---|---|
| default | 10 | 2.689e-13 | 5.389e-13 | 2.451e-13 | 2.806e-13 | 7.693e-14 | 3.080e-13 |
|  | 20 | 1.070e-12 | 7.111e-13 | 1.161e-12 | 6.813e-13 | 1.298e-12 | 6.731e-13 |
|  | 40 | 1.492e-12 | 2.950e-12 | 1.481e-12 | 2.955e-12 | 1.248e-12 | 1.988e-12 |
| `unequal` | 10 | 3.032e-13 | 4.854e-13 | 1.637e-13 | 2.584e-13 | 3.400e-13 | 6.903e-13 |
|  | 20 | 5.437e-13 | 4.913e-13 | 6.848e-13 | 6.368e-13 | 1.219e-12 | 1.169e-12 |
|  | 40 | 1.449e-12 | 2.903e-12 | 9.855e-13 | 2.929e-12 | 1.008e-12 | 2.998e-12 |
| bidomain | 10 | 1.660e-13 | 1.147e-13 | 1.470e-13 | 9.779e-14 | 1.206e-13 | 6.607e-14 |
|  | 20 | 3.695e-13 | 1.957e-13 | 2.470e-13 | 8.419e-14 | 2.652e-13 | 9.468e-14 |
|  | 40 | 5.191e-13 | 8.134e-13 | 1.842e-13 | 4.604e-13 | 4.148e-13 | 7.231e-13 |
| parallel (2 ranks) | 20 | 3.778e-12 | 3.076e-12 | 3.591e-12 | 1.566e-12 | 2.062e-12 | 2.298e-12 |

Tetrahedra:

| configuration | N | surface dynamic | surface static | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|---|---|
| default (corrected, anisotropic) | 10 | 6.033e-11 | 5.395e-11 | 4.646e-11 | 4.142e-11 | 3.756e-11 | 3.455e-11 |
|  | 14 | 1.943e-10 | 1.975e-10 | 1.736e-10 | 1.759e-10 | 2.484e-10 | 2.536e-10 |
|  | 20 | 1.073e-10 | 1.085e-10 | 8.177e-11 | 8.257e-11 | 1.245e-10 | 1.252e-10 |
|  | 28 | 1.677e-10 | 1.692e-10 | 1.129e-10 | 1.137e-10 | 1.102e-10 | 1.099e-10 |
| `unequal` | 10 | 6.098e-11 | 5.197e-11 | 4.822e-11 | 4.305e-11 | 3.212e-11 | 3.344e-11 |
|  | 14 | 5.846e-10 | 5.921e-10 | 3.540e-10 | 3.571e-10 | 2.532e-10 | 2.536e-10 |
|  | 20 | 2.842e-10 | 2.834e-10 | 3.782e-10 | 3.766e-10 | 2.704e-10 | 2.716e-10 |
|  | 28 | 1.520e-10 | 1.529e-10 | 1.006e-10 | 1.008e-10 | 1.105e-10 | 1.105e-10 |
| `isotropic` | 10 | 2.589e-12 | 3.040e-12 | 1.310e-12 | 1.831e-12 | 1.047e-12 | 1.672e-12 |
|  | 14 | 4.013e-12 | 3.864e-12 | 3.337e-12 | 4.007e-12 | 1.783e-12 | 1.935e-12 |
|  | 20 | 4.890e-12 | 3.646e-12 | 2.693e-12 | 2.819e-12 | 3.183e-12 | 3.372e-12 |
| bidomain | 10 | 1.610e-11 | 1.588e-11 | 1.351e-11 | 1.299e-11 | 1.098e-11 | 1.081e-11 |
|  | 14 | 1.041e-10 | 1.050e-10 | 6.000e-11 | 5.992e-11 | 1.020e-10 | 1.006e-10 |
|  | 20 | 5.207e-11 | 5.157e-11 | 4.039e-11 | 3.970e-11 | 4.418e-11 | 4.488e-11 |
|  | 28 | 5.037e-11 | 5.036e-11 | 3.682e-11 | 3.687e-11 | 3.943e-11 | 3.937e-11 |
| parallel (2 ranks) | 20 | 1.104e-10 | 1.075e-10 | 8.564e-11 | 8.465e-11 | 1.281e-10 | 1.303e-10 |

- **Hexahedra.** The gap on the body surface is at round-off and rises
  slightly with the cell count. It is of the same size for `unequal`, for
  the bidomain myocardium and on two ranks. Over all groups the gap lies
  between 6.6e-14 and 3.8e-12.
- **Tetrahedra.** The gap on the body surface is at the level of the
  fixed-point tolerances at every level and for every configuration, with
  no trend in `N`; over all groups it lies between 1.1e-11 and 5.9e-10.
  With only the non-orthogonal correction (`isotropic`) the fixed points
  take fewer passes (lead fields 32 to 59, direct solves 28 to 53, against
  40 to 77 and 39 to 84 for `default` and `unequal`) and the gap is
  1.0e-12 to 4.9e-12.
- **Parallel.** The parallel channels agree with the serial ones to 8.1e-12
  (absolute) on both meshes, at every electrode.
- **Probes.** The interior `torso` and `near` probes have gaps of the size
  of the surface group, `near` included, in torso cells that touch the
  heart. Electrodes inside the heart are not supported (case
  [README](../../../README.md)).

## Fixed-point tolerance control

`fixedPointTight_tet_N20` repeats the default tetrahedral run at `N = 20`
with both fixed points stopped at `1e-12` instead of `1e-10`:

| tolerance | surface dynamic | surface static | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|---|
| `1e-10` | 1.073e-10 | 1.085e-10 | 8.177e-11 | 8.257e-11 | 1.245e-10 | 1.252e-10 |
| `1e-12` | 4.030e-12 | 3.111e-12 | 2.137e-12 | 1.478e-12 | 2.208e-12 | 1.956e-12 |

Static Linf errors of the lead-field channel against the exact potential,
and of the direct static field over the surface cells (the direct channel
has the same electrode values):

| tolerance | LF surface | `phiStatic_surface` RMS | `phiStatic_surface` max | LF torso | LF near |
|---|---|---|---|---|---|
| `1e-10` | 1.5989e-03 | 1.0391e-03 | 1.6101e-03 | 1.6008e-03 | 1.6392e-03 |
| `1e-12` | 1.5989e-03 | 1.0391e-03 | 1.6101e-03 | 1.6008e-03 | 1.6392e-03 |

With the tighter tolerance the gap drops by more than an order of
magnitude in every group, to the round-off level of the hexahedra and of
`isotropic`, the errors against the exact potential are identical, and the
lead-field fixed points take about a third more passes. The tetrahedral gap
is set by the fixed-point tolerance, not by the transposed operator.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
