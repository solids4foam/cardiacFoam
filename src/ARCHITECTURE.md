# cardiacFoam Source Architecture

This document summarizes the libraries under `src/` as they exist in the
current tree. It is intentionally code-accurate rather than aspirational.

## Top-level layout

```text
src/
├── genericWriter/        Shared I/O and stimulus helpers
├── ionicModels/          Runtime-selectable cell-level ionic models
├── verificationModels/   Verification and manufactured-solution models
├── activeTensionModels/  Runtime-selectable active-tension ODE models
├── couplingModels/       Shared electromechanical signal interfaces
└── electroModels/        Spatial electrophysiology domains, solvers, couplers
```

Build order from `src/Allwmake`:

```text
genericWriter → ionicModels → verificationModels → activeTensionModels → electroModels
```

## Libraries

### `genericWriter` — `libgenericWriter`

Shared helpers for:

- ionic-model exports and trace writing
- ECG and Purkinje time-series output
- active-tension output support
- single-cell and monodomain stimulus parsing/evaluation

Main files:

- `ionicModelIO`
- `ionicVariableCompatibility`
- `stimulusIO`
- `activeTensionIO`
- `ecgModelIO`
- `purkinjeModelIO`

### `ionicModels` — `libionicModels`

Runtime-selectable cellular electrophysiology models. The base class is
`Foam::ionicModel`, with concrete models registered through
`addToRunTimeSelectionTable`.

This library is used by:

- `MyocardiumDomain` reaction-diffusion workflows
- `ConductionSystemDomain` graph-based workflows
- `singleCellSolver`

The `ionicModel/` subfolder now contains both the classic base/factory code and
batched or GPU-oriented support headers.

### `verificationModels` — `libverificationModels`

Verification infrastructure for spatial electrophysiology and ECG workflows.

Main layers:

- `electroVerification/` — base `electroVerificationModel`
- `monodomainVerification/`
- `bidomainVerification/`
- `ecgVerification/` — base `ecgVerificationModel` and ECG-specific verifiers

These are not ionic models. They are separate runtime-selected verifier
families that hook into myocardium or ECG workflows.

### `activeTensionModels` — `libactiveTensionModels`

Runtime-selectable active-tension models driven by an upstream
`ElectromechanicalSignalProvider`.

Current concrete models:

- `GoktepeKuhl`
- `NashPanfilov`

The implementation is integration-point based and scalar-state based. It is not
a tensor-mechanics framework on its own.

### `couplingModels`

Currently this folder is small and contains shared signal-side coupling
contracts, not the staged electro-domain couplers.

Current contents:

- `common/electromechanicalSignalProvider.H`

The staged Purkinje, ECG, and bath-style electro couplers live under
`src/electroModels/electroCouplers/`, not here.

### `electroModels` — `libelectroModels`

The main spatial electrophysiology stack. It contains:

- top-level orchestration in `core/`
- domain state owners in `electroDomains/`
- numerical solver kernels in `myocardiumModels/`, `conductionSystemModels/`,
  and `ecgModels/`
- staged inter-domain couplers in `electroCouplers/`

Current top-level electro entry is selected from `myocardiumSolver` in
`electroProperties`. The assembled multi-domain wrapper is
`electrophysiologyModel`, registered under:

- `monodomainSolver`
- `bidomainSolver`
- `eikonalSolver`

`singleCellSolver` is not a separate `src/` library. It is compiled inside
`electroModels/myocardiumModels/`.

## Reading guides

### Understand the current electro runtime path

1. `src/electroModels/README.md`
2. `src/electroModels/core/README.md`
3. `src/electroModels/core/ARCHITECTURE.md`

### Trace Purkinje-to-myocardium coupling

1. `src/electroModels/electroDomains/README.md`
2. `src/electroModels/conductionSystemModels/README.md`
3. `src/electroModels/electroCouplers/README.md`

### Add or debug an ionic model

1. `src/ionicModels/README.md`
2. `src/ionicModels/IONIC_MODEL_ARCHITECTURE.md`

### Understand verification hooks

1. `src/verificationModels/README.md`
2. `src/verificationModels/VERIFICATION_MODELS_ARCHITECTURE.md`

### Get the full electrophysiology picture

1. `src/electroModels/README.md`
2. `src/electroModels/ARCHITECTURE.md`
3. `src/electroModels/core/README.md`
4. `src/electroModels/core/ARCHITECTURE.md`
5. folder-level READMEs for the subsystem you are changing

