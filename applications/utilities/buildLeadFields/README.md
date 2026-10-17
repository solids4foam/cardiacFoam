# buildLeadFields

Builds the lead fields that the `leadFieldECG` ecgSolver records with. For
each electrode `e` of a `leadFieldECG` domain it solves the discrete adjoint
problem `A^T z_e = e_e`:

- `A` is the forward heart+torso operator of `passiveConductor`, with the
  conductivities, interface option, seal flag and schemes the case defines.
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
is the fixed point `M z = e_e - C^T z`.

The build checks, and stops with a fatal error above each tolerance:

| Check | Test | Tolerance |
|---|---|---|
| linearity | `C(r1 + 2 r2) = C r1 + 2 C r2` | 1e-11 |
| reproduction | probed `C` equals assembled `C` on another field | 1e-11 |
| transpose | `<z, C psi> = <C^T z, psi>` | 1e-12 |
| symmetry | `x^T M y = y^T M x`, coupled interfaces included | 1e-12 |
| residual | `e_e - M z - C^T z`, for each lead field | 1e-8 |
| pin | `z_r d_r = 1` at a pinned reference cell | 1e-8 |

Each value is written to the log (`leadFieldAdjoint: <check> <value>`). The
checks refuse nonlinear schemes (for example a `limited` snGrad or a limited
gradient). A cached `grad(phiE)` or `grad(phiE.component(0))` is fatal.

## Inputs

- `constant/[region/]electroProperties`: the `<myocardiumSolver>Coeffs`
  with the `leadFieldECG` entry of `ecgDomains`, its electrodes, Gi and Ge,
  and the torso conductor (the `bathPotentialDomain` block or the entry's
  `torso { ... }`).
- The mesh, the myocardium `cellZone`, the torso cellZones and the torso
  conductivity field (and the Gi/Ge fields with `conductivitySource field`).
- `system/fvSchemes`: the forward's own entries. The unknown is named
  `phiE` on the base mesh and on the heart sub-mesh, so the same
  `laplacian(...,phiE)` and `grad(phiE)` entries apply.
- `system/fvSolution`: `solvers { leadField { ... } leadFieldFinal { ... } }`,
  both required. The fixed point is solved until the relative change of `z`
  between passes is below 1e-10 (at least 3 passes, the last with
  `leadFieldFinal`; fatal after 100). The
  passes solve with `leadField` until the change first falls below 1e-10,
  then with `leadFieldFinal` until it falls below 1e-10 again.
  - `leadField` may use a loose `relTol` (e.g. 0.01).
  - `leadFieldFinal` must use `relTol 0` and a `tolerance` well below 1e-10.

  ```
  leadField
  {
      solver          PCG;
      preconditioner
      {
          preconditioner  GAMG;
          smoother        DICGaussSeidel;
      }
      tolerance       1e-15;
      relTol          0.01;
  }

  leadFieldFinal
  {
      $leadField;
      relTol          0;
  }
  ```

## Outputs

- `constant/[region/]leadField_<electrode>`: one `volScalarField` on the base
  mesh per electrode, in `1/(conductivity·length)`.
- `constant/[region/]leadFieldSignature`: the conductor's
  decomposition-independent signature (`operator adjoint`, cell counts,
  volume sums, conductivity sums, myocardium solver, interface option, seal
  flag, sealing patch types, bidomain heart trace, geometry scheme, the
  fvSchemes entries that resolve the schemes of `phiE` on the base and heart
  meshes, electrodes, ground and surface-current patches, reference).
  `leadFieldECG` compares it with the case at load.

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
| `-ecgDomain <name>` | the only `leadFieldECG` entry | `ecgDomains` entry to build; fatal when none or several exist without it |

## Example

```bash
buildLeadFields -case ./bathTorso
mpirun -np 2 buildLeadFields -parallel -ecgDomain leadField
```

Rebuild the lead fields whenever the mesh, conductivities, interface option,
seal flag, electrodes or schemes change.
