# verificationModels

This folder builds `libverificationModels`, the verification library for
electrophysiology and ECG workflows.

## Current contents

```text
src/verificationModels/
├── electroVerification/      # Base electro verifier family
├── monodomainVerification/   # Monodomain manufactured/reference verifiers
├── bidomainVerification/     # Bidomain manufactured/reference verifiers
├── ecgVerification/          # ECG verification family
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
  Base runtime-selection layer for myocardium-side verification hooks
- `ecgVerificationModel`
  Base runtime-selection layer for ECG-side verification hooks

## Concrete families

- `monodomainVerification/`
  Manufactured monodomain references and verifiers
- `bidomainVerification/`
  Manufactured bidomain references and verifiers
- `ecgVerification/`
  ECG verification helpers such as pseudo-ECG manufactured verification

## What this folder does not own

This folder does not own the myocardium or ECG solvers themselves. It provides
verification-side models that are called from those workflows.

