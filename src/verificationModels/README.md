# verificationModels

This folder builds `libverificationModels`, the verification library for
electrophysiology and ECG workflows.

## Current contents

```text
src/verificationModels/
├── monodomainVerification/   # Monodomain manufactured/reference verifiers
├── bidomainVerification/     # Bidomain manufactured/reference verifiers
├── bathBidomainVerification/ # Bath-bidomain manufactured/reference verifiers
├── eikonalVerification/      # Eikonal manufactured/reference verifiers
├── ecgVerification/          # ECG verifier family (concrete verifiers only)
├── electromechanicsVerification/ # Electromechanics verifiers
├── Make/
└── README.md
```

## Purpose

The code in this folder provides runtime-selected verification models and
reference helpers used to validate:

- monodomain workflows
- bidomain workflows
- ECG workflows
- selected single-cell manufactured cases

## Main abstractions

- `electroVerificationModel`
  Abstract base for myocardium-side verification hooks. Compiled in
  `electroModels/core/verificationModels/`.
- `ecgVerificationModel`
  Abstract base for ECG-side verification hooks. Compiled in
  `electroModels/core/verificationModels/`.
- `eikonalVerificationModel`
  Abstract base for eikonal activation-time verification hooks. Compiled in
  `electroModels/core/verificationModels/`.

## Concrete families

- `monodomainVerification/`
  Manufactured monodomain references and verifiers
- `bidomainVerification/`
  Manufactured bidomain references and verifiers
- `eikonalVerification/`
  Manufactured eikonal references and verifiers
- `ecgVerification/`
  ECG verification helpers such as pseudo-ECG manufactured verification

Registered verifier types include:

- `manufacturedFDAMonodomainVerifier`
- `manufacturedFDABidomainVerifier`
- `manufacturedFDABathBidomainVerifier`
- `singleCellManufacturedFDABidomainVerifier`
- `manufacturedEikonalVerifier`
- `pseudoECGManufacturedVerifier`
- `bathECGManufacturedVerifier`

## What this folder does not own

This folder does not own the myocardium or ECG solvers themselves. It provides
verification-side models that are called from those workflows.
