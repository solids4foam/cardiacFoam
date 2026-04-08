# ecgModels

This directory contains runtime-selectable ECG evaluation kernels used by
`ECGDomain`. These models are downstream consumers of myocardium state: they
read the finalized electrical solution after the tissue advance and produce ECG
signals or derived potentials without feeding current back into the tissue.

## Current contents

```text
src/electroModels/ecgModels/
├── pseudoECGSolver/   # Pseudo-ECG evaluation from myocardium state
└── README.md
```

## Available models

- `PseudoECGSolver`
  - Registered as `pseudoECG`.
  - Implements the current ECG evaluation kernel used by `ECGDomain`.
  - Reads upstream myocardium state through `ECGDomain`, which in turn depends
    on the `electroStateProvider` interface exposed by the primary domain.

## Execution role

`ECGDomain` is a downstream domain in the `electrophysicsSystem`:

- it advances after the myocardium
- it consumes already-updated tissue state
- it does not inject source terms back into the myocardium

See [../electroDomains/README.md](../electroDomains/README.md) for the
domain-level contract and [../core/ELECTROMODEL_ORCHESTRATION.md](../core/ELECTROMODEL_ORCHESTRATION.md)
for the timestep sequence.
