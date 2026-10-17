# buildLeadFields

Builds the lead fields that the `leadFieldECG` ecgSolver records with. For
each electrode `e` of a `leadFieldECG` domain it solves the reciprocal
problem `A z_e = δ_e`:

- `A` is the forward heart+torso operator of the bath bidomain
  (`passiveConductor`), with the conductivities, interface option and seal
  flag the case defines.
- `δ_e` is a unit current in the cell that contains the electrode
  (`1/V` in that cell).
- Ground patches are `fixedValue 0`; every other non-coupled patch is
  `zeroGradient`. Without ground patches, the `phiERefPoint` cell is pinned
  to 0.

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
  both required. The operator is reassembled and solved until the relative
  change of `z` between passes is below 1e-10 (at least 2 passes; fatal after
  100). The passes solve with `leadField` until the change first falls below
  1e-10, then with `leadFieldFinal` until it falls below 1e-10 again.
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
  decomposition-independent signature (cell counts, volume sums, conductivity
  sums, interface option, seal flag, sealing patch types, electrodes, ground
  and surface-current patches, reference). `leadFieldECG` compares it with
  the case at load.

In parallel (`-parallel`), each processor writes its own
`processorN/constant/[region/]leadField_<electrode>`. A parallel ECG run
therefore needs lead fields built with the same decomposition (run
`buildLeadFields -parallel` after `decomposePar`) or decomposed from a serial
build.

## Options

| Option | Default | Description |
|---|---|---|
| `-region <name>` | `region0` | Mesh region |
| `-ecgDomain <name>` | the only `leadFieldECG` entry | `ecgDomains` entry to build; fatal when none or several exist without it |

## Example

    buildLeadFields -case ./bathTorso
    mpirun -np 2 buildLeadFields -parallel -ecgDomain leadField

Rebuild the lead fields whenever the mesh, conductivities, interface option,
seal flag, electrodes or schemes change.
