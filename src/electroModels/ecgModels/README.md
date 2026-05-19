# ecgModels

This directory contains runtime-selectable ECG evaluation kernels used by
`ecgDomain`. These models are downstream consumers of myocardium state: they
read the finalized electrical solution after the tissue advance and produce ECG
signals or derived potentials without feeding current back into the tissue.

## Current contents

```text
src/electroModels/ecgModels/
├── pseudoECGSolver/
│   ├── pseudoECGSolver.H
│   └── pseudoECGSolver.C
├── torsoECG/
│   ├── torsoECG.H
│   └── torsoECG.C
└── README.md
```

## Available models

**Concrete solver implementations:**

- **`pseudoECGSolver`** (ECG post-processor)
  - Registered as `pseudoECG`.
  - Computes pseudo-ECG signals using the Gima-Rudy dipole model.
  - Reads upstream myocardium state through `ecgDomain`.
  - Abstract interface: `electroDomains/ecgDomain/ecgSolver.H/C`

- **`torsoECG`** (electrode sampler on the unified bath potential)
  - Registered as `torsoECG`.
  - Samples the globally solved `phiE` field at electrode positions on the
    union (heart + bath) mesh; parallel-safe via list reduction.
  - The global `phiE` solve is owned by `extracellularPotentialDomain` — see
    [../electroDomains/extracellularPotentialDomain/](../electroDomains/extracellularPotentialDomain/).
  - Selected by routing the ECG-domain state provider to the configured
    `bathPotentialDomain` inside `bidomainSolverCoeffs` (done in
    `electrophysicsSystemBuilder::configureECGDomains`).
  - Abstract interface: `electroDomains/ecgDomain/ecgSolver.H/C`

## Architectural pattern

- **Abstract solver interface** lives in the domain folder:
  - `electroDomains/ecgDomain/ecgSolver.H/C`

- **Concrete solver implementations** live here in `ecgModels/`:
  - `pseudoECGSolver/`
  - `torsoECG/`

## Execution role

`ecgDomain` is a downstream domain in the `electrophysicsSystem`:

- it advances after the myocardium
- it consumes already-updated tissue state
- it does not inject source terms back into the myocardium

See [../electroDomains/README.md](../electroDomains/README.md) for the
domain-level contract and [../core/ARCHITECTURE.md](../core/ARCHITECTURE.md)
for the timestep sequence.

`torsoECG` is part of the active orchestration path. It is wired by
`electrophysicsSystemBuilder` when an ECG domain selects `torsoECG` and a
`bidomainSolverCoeffs.bathPotentialDomain` block provides the global `phiE`
state.
