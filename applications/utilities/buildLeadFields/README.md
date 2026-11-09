# buildLeadFields

Builds the lead fields that `ecgSolver torsoECG` records with under
`torsoCoupling leadField`. For each electrode `e` of the `torsoECG` domain it
solves the discrete adjoint
problem `A^T z_e = e_e`:

- `A` is the forward heart+torso operator of `passiveConductor`, with the
  conductivities, interface option, seal flag and schemes the case defines
  (with `conormalHarmonic`, `nInterfaceGradientCorrectors` of the `phiE`
  entry of `solvers` in `system/fvSolution`).
- `e_e` is the unit vector of the cell that contains the electrode.
- Ground patches are `fixedValue 0`; every other non-coupled patch is
  `zeroGradient`. Without ground patches, the `phiERefPoint` cell is pinned
  to 0.

Then `z_e^T b` equals the electrode value of the converged forward solve
`A phi = b` of the same conductor on any mesh.

`A = M + C`: `M` is the implicit matrix, assembled once and symmetric; `C` is
the explicit part (cross-diffusion, non-orthogonal and interface
corrections). `C` is stored as a sparse matrix by probing the same assembly
with one field per colour of a distance-2 colouring of its stencil
(`probedOperator`), so `C^T` is exact for any linear scheme. Each lead field
solves `(M + C^T) z = e_e` with the assembled matrix through PETSc
(`conductorLinearSolve`, solids4foam's `foamPetscSnesHelper`); the matrix
and its preconditioner are set up once for all electrodes. A PETSc solve that
does not converge is fatal.

The build checks, and stops with a fatal error above each tolerance:

| Check | Test | Tolerance |
|---|---|---|
| linearity | `C(r1 + 2 r2) = C r1 + 2 C r2` | 1e-11 |
| reproduction | probed `C` equals assembled `C` on another field | 1e-11 |
| transpose | `<z, C psi> = <C^T z, psi>` | 1e-12 |
| symmetry | `x^T M y = y^T M x`, coupled interfaces included | 1e-12 |
| residual | `e_e - M z - C^T z`, for each lead field | `residualTolerance` |
| pin | `z_r d_r = 1` at a pinned reference cell | `residualTolerance` |

The first four are round-off checks of exactly linear operators. The residual
and the pin are the true residual of the solve, checked against the
`residualTolerance` of the `leadField` entry (below). Each value is
written to the log (`leadFieldAdjoint: <check> <value>`). The checks refuse
nonlinear schemes (for example a `limited` snGrad or a limited gradient). A
cached `grad(phiE)` or `grad(phiE.component(0))` is fatal.

## Inputs

- `constant/[region/]electroProperties`: `torsoCoupling leadField`, the
  top-level `torso` block (the torso conductor), the `torsoECG` entry of the
  top-level `ecgDomains` block (its electrodes, and Gi and Ge for a
  `monodomainSolver` heart) and the `<myocardiumSolver>Coeffs` (the heart,
  and Gi and Ge for a `bidomainSolver` heart). Its keywords are checked as
  `cardiacFoam` checks them.
- The mesh, the myocardium `cellZone`, the torso cellZones and the torso
  conductivity field (and the Gi/Ge fields with `conductivitySource field`).
- `system/fvSchemes`: the forward's own entries. The unknown is named
  `phiE` on the base mesh and on the heart sub-mesh, so the same
  `laplacian(...,phiE)` and `grad(phiE)` entries apply.
- `system/fvSolution`: the solver entry `leadField` with `solver petsc`,
  `residualTolerance` (required) and the PETSc `options`. The build needs
  cardiacFoam built with PETSc and the full solids4foam; PETSc runs with
  floating-point traps masked. Cyclic patches are not supported. The
  lead-field ECG verifiers (`manufacturedLeadFieldECGVerifier`,
  `manufacturedShellLeadFieldECGVerifier`) solve their direct forward
  problem the same way with the entry `leadFieldDirect`.

  ```
  leadField
  {
      solver              petsc;
      residualTolerance   1e-10;
      options
      {
          snes_type           ksponly;
          ksp_type            gmres;
          ksp_gmres_restart   "200";
          ksp_pc_side         right;
          ksp_norm_type       unpreconditioned;
          ksp_rtol            "1e-12";
          ksp_atol            "0";
          ksp_max_it          "5000";
          pc_type             gamg;
      }
  }
  ```

## Outputs

- `constant/[region/]leadField_<electrode>`: one `volScalarField` on the base
  mesh per electrode, in `1/(conductivity·length)`.
- `constant/[region/]leadFieldSignature`: the conductor's
  decomposition-independent signature (`operator adjoint`, cell counts,
  volume sums, conductivity sums, myocardium solver, interface option (and
  `nInterfaceGradientCorrectors` with `conormalHarmonic`), seal
  flag, sealing patch types, bidomain heart trace, geometry scheme, the
  fvSchemes entries that resolve the schemes of `phiE` on the base and heart
  meshes, electrodes, ground and surface-current patches, reference).
  `torsoECG` compares it with the case at load.
- `constant/[region/]leadFieldBuildRecord`: how the lead fields were built:
  the `leadField` solver entry, and for each electrode the PETSc
  iterations, the converged reason and a fingerprint, the sums of
  `z V f` over the cells for `f = 1, x, y, z, |x|²`. The signature holds
  totals, which a renumbered or redistributed mesh keeps; the fingerprint
  pairs each value with its cell. `torsoECG` checks each lead field's
  fingerprint at load and logs its iterations and converged reason.

In parallel (`-parallel`), each processor writes its own
`processorN/constant/[region/]leadField_<electrode>`. A parallel ECG run
therefore needs lead fields built with the same decomposition (run
`buildLeadFields -parallel` after `decomposePar`), or a serial build
decomposed with `decomposePar -constant` (the lead fields are in
`constant`).

## Options

| Option | Default | Description |
|---|---|---|
| `-region <name>` | `region0` | Mesh region |
| `-ecgDomain <name>` | the only `torsoECG` entry | `ecgDomains` entry to build; fatal when none or several exist without it |

## Example

```bash
buildLeadFields -case ./bathTorso
mpirun -np 2 buildLeadFields -parallel -ecgDomain leadField
```

Rebuild the lead fields whenever the mesh, conductivities, interface option,
interface-gradient passes, seal flag, electrodes or schemes change.
