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

All runs use the bidomain myocardium (`-myocardium bidomain`). The default
configuration is the corrected, anisotropic case of
[leadFieldConvergence](../leadFieldConvergence/README.md): hexahedra
`N = 10, 20, 40`, tetrahedra `N = 10, 14, 20, 28`; that study also runs
`unequal` (Gi and Ge with different anisotropy ratios) on the same levels.
This study adds:

| configuration | options | levels |
|---|---|---|
| `isotropic_tet` | `-mesh tet -conductivity isotropic` (non-orthogonal correction only) | 10, 14, 20 |
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
    bash results/$name/Allrun "$@" -myocardium bidomain -endTime 0.02
}
for n in 10 14 20
do
    run isotropic_tet_N$n -mesh tet -N $n -conductivity isotropic
done
run parallel_hex_N20 parallel -mesh hex -N 20
run parallel_tet_N20 parallel -mesh tet -N 20
run solveTight_tet_N20 -mesh tet -N 20 -solveTolerance 1e-14
python3 ../summarize.py ../leadFieldConvergence/results results
```

The omnidriver specs for this study are `sweep_hex_reciprocity.json` and
`sweep_tet_reciprocity.json`, for the record `manufacturedLeadFieldECG`, which
omnidriver does not have. They name the `Allrun` options as axes
(`conductivity`, `mesh`, `parallel`, with `myocardium bidomain` in the base; `-N` is `numberCells` on hexahedra and `tetNumberCells` on tetrahedra),
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
| default (anisotropic) | 10 | 6.853e-13 | 6.226e-13 | 3.758e-13 | 3.136e-13 | 4.973e-13 | 4.304e-13 |
|  | 20 | 3.280e-13 | 3.761e-13 | 4.179e-13 | 4.651e-13 | 3.828e-13 | 4.322e-13 |
|  | 40 | 9.181e-13 | 1.078e-12 | 9.435e-13 | 1.104e-12 | 7.422e-13 | 9.060e-13 |
| `unequal` | 10 | 5.852e-13 | 5.368e-13 | 5.782e-13 | 5.298e-13 | 2.282e-13 | 1.792e-13 |
|  | 20 | 8.284e-13 | 6.912e-13 | 6.595e-13 | 7.148e-13 | 5.781e-13 | 6.315e-13 |
|  | 40 | 7.449e-13 | 6.041e-13 | 3.921e-13 | 2.495e-13 | 5.824e-13 | 4.327e-13 |
| parallel (2 ranks) | 20 | 1.343e-12 | 8.325e-13 | 1.223e-12 | 1.836e-12 | 8.678e-13 | 8.392e-13 |

Tetrahedra:

| configuration | N | surface dynamic | surface static | torso dynamic | torso static | near dynamic | near static |
|---|---|---|---|---|---|---|---|
| default (anisotropic) | 10 | 4.689e-13 | 4.647e-13 | 1.566e-13 | 1.432e-13 | 3.299e-13 | 3.305e-13 |
|  | 14 | 7.558e-13 | 6.960e-13 | 7.136e-13 | 6.653e-13 | 2.546e-13 | 1.999e-13 |
|  | 20 | 8.129e-13 | 1.014e-12 | 7.295e-13 | 6.749e-13 | 1.259e-12 | 1.479e-12 |
|  | 28 | 7.171e-13 | 1.093e-12 | 4.386e-13 | 6.517e-13 | 3.033e-13 | 6.904e-13 |
| `unequal` | 10 | 3.139e-13 | 3.295e-13 | 3.207e-13 | 1.872e-13 | 2.962e-13 | 1.688e-13 |
|  | 14 | 6.891e-13 | 7.764e-13 | 5.513e-13 | 6.403e-13 | 3.675e-13 | 4.443e-13 |
|  | 20 | 1.117e-12 | 1.129e-12 | 7.362e-13 | 7.411e-13 | 6.963e-13 | 7.182e-13 |
|  | 28 | 7.932e-13 | 1.041e-12 | 1.142e-12 | 8.665e-13 | 7.835e-13 | 1.220e-12 |
| `isotropic` | 10 | 8.563e-13 | 9.653e-13 | 8.839e-13 | 9.893e-13 | 3.864e-13 | 4.990e-13 |
|  | 14 | 1.063e-12 | 9.433e-13 | 2.259e-12 | 2.128e-12 | 4.784e-13 | 3.859e-13 |
|  | 20 | 2.860e-12 | 2.719e-12 | 9.756e-13 | 8.214e-13 | 9.446e-13 | 7.757e-13 |
| parallel (2 ranks) | 20 | 1.881e-12 | 1.809e-12 | 7.520e-13 | 7.490e-13 | 1.534e-12 | 1.165e-12 |

- **Hexahedra.** The gap is at round-off in every group, for the default
  and `unequal` conductivities and on two ranks: between 1.8e-13 and
  1.9e-12.
- **Tetrahedra.** The gap is at round-off as well, with no trend in `N` and
  no dependence on the configuration: between 1.4e-13 and 1.9e-12, and
  3.9e-13 to 2.9e-12 with only the non-orthogonal correction
  (`isotropic`). The lead fields take 25 to 42 GMRES iterations and the
  direct solves 27 to 45.
- **Parallel.** The parallel channels agree with the serial ones to 3.1e-12
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
| `1e-12` | 8.129e-13 | 1.014e-12 | 7.295e-13 | 6.749e-13 | 1.259e-12 | 1.479e-12 |
| `1e-14` | 4.637e-13 | 3.750e-13 | 1.114e-13 | 2.716e-13 | 1.877e-13 | 3.249e-13 |

Static Linf errors of the lead-field channel against the exact potential,
and of the direct static field over the surface cells (the direct channel
has the same electrode values):

| `ksp_rtol` | LF surface | `phiStatic_surface` RMS | `phiStatic_surface` max | LF torso | LF near |
|---|---|---|---|---|---|
| `1e-12` | 1.5997e-03 | 1.0403e-03 | 1.6104e-03 | 1.6014e-03 | 1.6311e-03 |
| `1e-14` | 1.5997e-03 | 1.0403e-03 | 1.6104e-03 | 1.6014e-03 | 1.6311e-03 |

With the tighter tolerance the gap drops by a factor of 2 to 7, the errors
against the exact potential are identical to every printed digit, and the
lead fields take 38 to 40 iterations instead of 32 to 34. The solve
tolerance does not reach the reported errors.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
