# Verification Models Architecture

This document describes the verification architecture that exists in the
current tree. It focuses on the runtime-selected verifier layers and the
manufactured/reference helpers they use.

## Overview

`src/verificationModels` is split into two verifier families plus the shared
manufactured/reference code they use:

- `electroVerificationModel`
  abstract base class — defined in `electroModels/core/verificationModels/`.
- `ecgVerificationModel`
  abstract base class — defined in `electroModels/core/verificationModels/`.
- `eikonalVerificationModel`
  abstract base class — defined in `electroModels/core/verificationModels/`.

## Directory layout

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

## Myocardium-side verification

### `electroVerificationModel`

Defined in:

- `electroModels/core/verificationModels/electroVerificationModel.H` ← **authoritative location**
- `electroModels/core/verificationModels/electroVerificationModel.C`

Role:

- runtime-selection base for myocardium-side verifiers
- constructed from dictionary input
- provides verification hooks that spatial workflows can call before or after
  the main solve

Current concrete verifiers include:

- `manufacturedFDAMonodomainVerifier`
- `manufacturedFDABidomainVerifier`
- `manufacturedFDABathBidomainVerifier`

### Manufactured/reference helpers

Shared analytical helpers live beside the concrete verifiers:

- `monodomainVerification/manufacturedFDAReference.H`
- `bidomainVerification/manufacturedFDABidomainReference.H`

These provide exact/reference fields and helper formulas used by the concrete
verification models.

## Eikonal verification

### `eikonalVerificationModel`

Defined in:

- `electroModels/core/verificationModels/eikonalVerificationModel.H` ← **authoritative location**
- `electroModels/core/verificationModels/eikonalVerificationModel.C`

Role:

- runtime-selection base for eikonal activation-time verifiers
- constructed from dictionary input
- provides verification hooks that eikonal spatial workflows can call before or after
  the main solve

Current concrete verifier:

- `manufacturedEikonalVerifier`

### Eikonal manufactured/reference helpers

Shared analytical helpers live beside the concrete verifiers:

- `eikonalVerification/manufacturedEikonalReference.H`

These provide exact/reference fields and helper formulas used by the concrete
verification models.

## ECG-side verification

### `ecgVerificationModel`

Defined in:

- `electroModels/core/verificationModels/ecgVerificationModel.H` ← **authoritative location**
- `electroModels/core/verificationModels/ecgVerificationModel.C`

Role:

- runtime-selection base for ECG verification
- built around an upstream `electroStateProvider`
- owns ECG-specific verification requirements and lifecycle

Current concrete verifier:

- `pseudoECGManufacturedVerifier`

## What this layer validates

The current verification code is used to validate:

- monodomain spatial workflows
- bidomain spatial workflows
- ECG post-processing workflows
- selected single-cell manufactured workflows

It does that through:

- manufactured/reference formulas
- error reporting and summary hooks
- runtime-selected verifier objects integrated into the main electro stack

This folder does not own the myocardium or ECG solvers. It provides verification-side models called from those workflows.
