# leadFieldReciprocity - leadFieldECG

Does the lead-field channel equal the direct forward solve of the same
conductor on the same `Vm`?

The direct solve gives `e^T A^-1 b`, sampled in the electrode cell, with `A`
the forward operator (implicit part `M` and explicit face terms `C`). The
lead field is the discrete adjoint
`z_e = A^-T e`, so the channel `z_e^T b = e^T A^-1 b` equals it for any `A`.
The adjoint matters because `A` is not symmetric on general meshes: the
non-orthogonal correction of the `corrected` snGrad and OpenFOAM's tensor
cross-diffusion `(Sf.G - (Sf.G.n) n) . grad` are explicit terms whose
transpose differs from them, and they vanish only on orthogonal faces whose
normal is an eigenvector of the conductivity (the Cartesian hexahedra with
diagonal conductivities here). `buildLeadFields` transposes them exactly
(probed explicit operator, see its README). Both solves use the assembled
matrix `M + C` (or `M + C^T`) through PETSc, GMRES preconditioned by GAMG to
`ksp_rtol 1e-12`, and stop on the true residual of the assembled system
(`residualTolerance 1e-10` of `leadField` and `leadFieldDirect` in
`system/fvSolution`); the solve tolerance sets the size of the gap.

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
| `solveTight_tet` | `-mesh tet -solveTolerance 1e-14` (`ksp_rtol`) | 20 |

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
run solveTight_tet_N20 -mesh tet -N 20 -solveTolerance 1e-14
python3 ../summarize.py ../leadFieldConvergence/results results
```

The omnidriver specs for this study are `sweep_hex_reciprocity.json` and
`sweep_tet_reciprocity.json`, for the record `manufacturedLeadFieldECG`, which
omnidriver does not have. They name the `Allrun` options as axes
(`myocardium`, `conductivity`, `mesh`, `parallel`; `-N` is `numberCells` on hexahedra and `tetNumberCells` on tetrahedra),
`-endTime` as the keys `system/controlDict:endTime` and `writeInterval`, and
`-solveTolerance` as the keys
`system/fvSolution:solvers.leadField.options.ksp_rtol` and
`solvers.leadFieldDirect.options.ksp_rtol`.

## Measured behaviour

Gap relative to the largest exact potential (1.43), as the largest Linf over
the samples among the electrodes of the group: `dynamic` is
`leadField-direct` on `Vm_h`, `static` is `leadFieldStatic-directStatic` on
the exact `Vm`. Tetrahedra: 14792, 39131, 111482 and 299957 cells at
`N = 10, 14, 20, 28`.

Hexahedra (corrected; the snGrad scheme has no effect there):

| configuration | N | surface dynamic | surface static | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|---|---|
| default (corrected, anisotropic) | 10 | 8.143e-13 | 5.632e-13 | 1.056e-12 | 7.809e-13 | 3.540e-13 | 5.771e-13 |
|  | 14 | 1.189e-12 | 1.099e-12 | 7.390e-13 | 6.864e-13 | 1.001e-12 | 8.683e-13 |
|  | 20 | 9.630e-13 | 1.367e-12 | 5.779e-13 | 1.003e-12 | 5.630e-13 | 9.887e-13 |
|  | 28 | 3.348e-12 | 2.944e-12 | 2.567e-12 | 2.221e-12 | 3.332e-12 | 2.776e-12 |
| `unequal` | 10 | 2.640e-13 | 3.086e-13 | 3.541e-13 | 3.362e-13 | 6.410e-13 | 5.811e-13 |
|  | 14 | 1.753e-12 | 6.842e-13 | 8.557e-13 | 7.150e-13 | 1.344e-12 | 3.746e-13 |
|  | 20 | 1.020e-12 | 1.265e-12 | 1.146e-12 | 1.054e-12 | 6.606e-13 | 9.684e-13 |
|  | 28 | 2.645e-12 | 2.425e-12 | 2.713e-12 | 3.158e-12 | 2.318e-12 | 2.691e-12 |
| `isotropic` | 10 | 1.479e-12 | 1.292e-12 | 3.101e-13 | 4.532e-13 | 1.326e-12 | 1.170e-12 |
|  | 14 | 1.661e-12 | 1.703e-12 | 1.438e-12 | 1.489e-12 | 1.216e-12 | 1.243e-12 |
|  | 20 | 3.837e-12 | 3.387e-12 | 8.827e-13 | 8.414e-13 | 1.221e-12 | 1.187e-12 |
| bidomain | 10 | 4.689e-13 | 4.647e-13 | 1.566e-13 | 1.432e-13 | 3.299e-13 | 3.305e-13 |
|  | 14 | 7.558e-13 | 6.960e-13 | 7.136e-13 | 6.653e-13 | 2.546e-13 | 1.999e-13 |
|  | 20 | 8.129e-13 | 1.014e-12 | 7.295e-13 | 6.749e-13 | 1.259e-12 | 1.479e-12 |
|  | 28 | 7.171e-13 | 1.093e-12 | 4.386e-13 | 6.517e-13 | 3.033e-13 | 6.904e-13 |
| parallel (2 ranks) | 20 | 3.925e-12 | 4.136e-12 | 4.906e-12 | 5.121e-12 | 2.488e-12 | 1.564e-12 |

- **Hexahedra.** The gap is at round-off in every group, for `unequal`, for
  the bidomain myocardium and on two ranks: between 1.9e-13 and 4.2e-12.
- **Tetrahedra.** The gap is at round-off as well, with no trend in `N` and
  no dependence on the configuration: between 1.4e-13 and 5.1e-12, and
  3.1e-13 to 3.8e-12 with only the non-orthogonal correction
  (`isotropic`). The lead fields take 25 to 44 GMRES iterations and the
  direct solves 27 to 47.
- **Parallel.** The parallel channels agree with the serial ones to 5.9e-12
  (absolute) on both meshes, at every electrode.
- **Probes.** The interior `torso` and `near` probes have gaps of the size
  of the surface group, `near` included, in torso cells that touch the
  heart. Electrodes inside the heart are not supported (case
  [README](../../../README.md)).

## Solve tolerance control

`solveTight_tet_N20` repeats the default tetrahedral run at `N = 20` with
both solves at `ksp_rtol 1e-14` instead of `1e-12`:

| `ksp_rtol` | surface dynamic | surface static | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|---|
| `1e-12` | 9.630e-13 | 1.367e-12 | 5.779e-13 | 1.003e-12 | 5.630e-13 | 9.887e-13 |
| `1e-14` | 4.860e-13 | 7.386e-13 | 4.249e-13 | 2.493e-13 | 4.564e-13 | 3.220e-13 |

Static Linf errors of the lead-field channel against the exact potential,
and of the direct static field over the surface cells (the direct channel
has the same electrode values):

| `ksp_rtol` | LF surface | `phiStatic_surface` RMS | `phiStatic_surface` max | LF torso | LF near |
|---|---|---|---|---|---|
| `1e-12` | 1.5989e-03 | 1.0391e-03 | 1.6101e-03 | 1.6008e-03 | 1.6392e-03 |
| `1e-14` | 1.5989e-03 | 1.0391e-03 | 1.6101e-03 | 1.6008e-03 | 1.6392e-03 |

With the tighter tolerance the gap drops by a factor of 2 to 4, the errors
against the exact potential are identical to every printed digit, and the
lead fields take 39 to 41 iterations instead of 33 to 36. The solve
tolerance does not reach the reported errors.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
