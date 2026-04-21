# Verification Models Architecture

This document describes the verification architecture that exists in the
current tree. It focuses on the runtime-selected verifier layers and the
manufactured/reference helpers they use.

## Overview

`src/verificationModels` is split into two verifier families plus the shared
manufactured/reference code they use:

- `electroVerificationModel`
  for myocardium-side verification hooks
- `ecgVerificationModel`
  for ECG-side verification hooks

The code here is not structured as `ionicModel` subclasses. That was an older
shape. The current implementation uses dedicated verifier base classes.

## Directory layout

```text
src/verificationModels/
├── electroVerification/      # Base electro verifier family
├── monodomainVerification/   # Monodomain manufactured/reference verifiers
├── bidomainVerification/     # Bidomain manufactured/reference verifiers
├── ecgVerification/          # ECG verifier family
├── Make/
└── README.md
```

## Myocardium-side verification

### `electroVerificationModel`

Defined in:

- `electroVerification/electroVerificationModel.H`
- `electroVerification/electroVerificationModel.C`

Role:

- runtime-selection base for myocardium-side verifiers
- constructed from dictionary input
- provides verification hooks that spatial workflows can call before or after
  the main solve

Current concrete verifiers include:

- `manufacturedFDAMonodomainVerifier`
- `manufacturedFDABidomainVerifier`
- `singleCellManufacturedFDABidomainVerifier`

### Manufactured/reference helpers

Shared analytical helpers live beside the concrete verifiers:

- `monodomainVerification/manufacturedFDAReference.H`
- `bidomainVerification/manufacturedFDABidomainReference.H`

These provide exact/reference fields and helper formulas used by the concrete
verification models.

## ECG-side verification

### `ecgVerificationModel`

Defined in:

- `ecgVerification/ecgVerificationModel.H`
- `ecgVerification/ecgVerificationModel.C`

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

## What this layer does not own

This folder does not own:

- the ionic model hierarchy
- the myocardium solver kernels
- the ECG solver kernels

It provides verification-side models and reference utilities that are called
from those workflows.
