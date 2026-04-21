# ecgModels

This directory contains runtime-selectable ECG evaluation kernels used by
`ECGDomain`. These models are downstream consumers of myocardium state: they
read the finalized electrical solution after the tissue advance and produce ECG
signals or derived potentials without feeding current back into the tissue.

## Current contents

```text
src/electroModels/ecgModels/
├── pseudoECGSolver/
│   ├── pseudoECGSolver.H
│   └── pseudoECGSolver.C
├── bidomainBathECGSolver/
│   ├── bidomainBathECGSolver.H
│   └── bidomainBathECGSolver.C
└── README.md
```

## Available models

**Concrete solver implementations:**

- **`PseudoECGSolver`** (ECG post-processor)
  - Registered as `pseudoECG`.
  - Computes pseudo-ECG signals using the Gima-Rudy dipole model.
  - Reads upstream myocardium state through `ECGDomain`.
  - Abstract interface: `electroDomains/ecgDomain/ecgSolver.H/C`

- **`BidomainBathECGSolver`** (bath extracellular potential solver)
  - Registered as `bidomainBathECG`.
  - Solves steady-state Laplacian: `∇·(σ_bath·∇φE) = -I_interface`
  - Reads myocardium transmembrane current through `BathDomain`.
  - Abstract interface: `electroDomains/bathDomain/bathECGSolver.H/C`

## Architectural pattern

- **Abstract solver interfaces** live in domain folders (`electroDomains/`):
  - `ecgDomain/ecgSolver.H/C`
  - `bathDomain/bathECGSolver.H/C`
  
- **Concrete solver implementations** live here in `ecgModels/`:
  - `pseudoECGSolver/`
  - `bidomainBathECGSolver/`

## Execution role

`ECGDomain` is a downstream domain in the `electrophysicsSystem`:

- it advances after the myocardium
- it consumes already-updated tissue state
- it does not inject source terms back into the myocardium

See [../electroDomains/README.md](../electroDomains/README.md) for the
domain-level contract and [../core/ARCHITECTURE.md](../core/ARCHITECTURE.md)
for the timestep sequence.

Bath-related solver code is still present in this folder, but bath is not part
of the active `core` orchestration path at the moment.
