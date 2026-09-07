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
├── electroModels/        Spatial electrophysiology domains, solvers, couplers
└── electroMechanicalModels/ Electromechanics wrappers built in full solids4foam mode

```

Build order from `src/Allwmake`:

```text
couplingModels lnInclude → genericWriter → ionicModels →
activeTensionModels → electroModels → verificationModels →
electroMechanicalModels (full mode only)

```

`etc/resolveSolids4Foam.sh` selects a built solids4foam installation when one
is available; otherwise it builds and uses the repository's lightweight
`physicsModel` compatibility library. Both modes build the electrophysiology
libraries through `verificationModels`. Only full mode builds
`electroMechanicalModels`; solid mechanics and coupled electromechanics are
therefore unavailable in lightweight mode.

## Runtime selection layers

| Layer | Public dictionary surface | Runtime owner |
|---|---|---|
| Top-level physics | `constant/physicsProperties`: `type electroModel` | `physicsModel::New()` constructs `electroModel` |
| Electro workflow | `constant/electroProperties`: `myocardiumSolver` | `electroModel::New()` selects a direct model or spatial wrapper |
| Spatial assembly | `<myocardiumSolver>Coeffs` | `electrophysiologyModel` and `electrophysicsSystemBuilder` assemble domains and couplers |
| Myocardium kernel | `monodomainSolver` or `bidomainSolver` | `myocardiumDomain` owns a runtime-selected `myocardiumSolver`; eikonal uses `eikonalMyocardiumDomain` |

The spatial names `monodomainSolver`, `bidomainSolver`, and `eikonalSolver` are
aliases for `electrophysiologyModel` in the parent `electroModel` table.
`singleCellSolver` registers directly in that parent table and does not create
the multi-domain system.

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

- `myocardiumDomain` reaction-diffusion workflows

- `conductionSystemDomain` graph-based workflows

- `singleCellSolver`

The `ionicModel/` subfolder now contains both the classic base/factory code and
batched or GPU-oriented support headers.

**Tissue heterogeneity support:** Ionic models can optionally configure spatial
heterogeneity of cellular phenotypes (endocardial, mid-myocardial, epicardial, or
open-ended named/scar regions) through the `ionicHeterogeneity` dictionary block,
in one of three modes — `transmuralBands` (fixed 3-zone), `namedRegions` (open,
field-range-based), or `cellZoneRegions` (mesh-topology-based) — plus an optional
`apexBaseBands` scaling composable with any mode. All 12 scalar models (via
`configuredIonicModel`) support this; batched models support it only when their
`supportedTissueTypes()` includes all three anatomical types (currently
`BuenoOrovioBatched`, `TNNPBatched`, `TWorldBatched`, `ToRORd_dynClBatched`). See
`src/ionicModels/README.md` for configuration details.

### `verificationModels` — `libverificationModels`

Verification infrastructure for spatial electrophysiology and ECG workflows.

This library depends on `electroModels` (compiled after it) because concrete
verifiers inherit from base classes (e.g., `electroVerificationModel`, `ecgVerificationModel`, `eikonalVerificationModel`, `couplingVerificationModel`, and `graphVerificationModel`)
which are defined in `electroModels/core/verificationModels/`.

Main layers:

- `monodomainVerification/`
- `bidomainVerification/`
- `ecgVerification/`
- `coupledVerification/`

### `activeTensionModels` — `libactiveTensionModels`

Runtime-selectable active-tension models driven by an upstream
`ElectromechanicalSignalProvider`.

Current concrete models:

- `NashPanfilov`

- `LandNiederer` (the original seven-state intact-human model)

- `LandNiedererTWorld`

`LandNiedererTWorld` ships the matching GPU-batched variant. The original
`LandNiederer` passive branch is diagnostic only because solid mechanics owns
the passive constitutive response.

### `couplingModels`

Currently this folder is small and contains shared signal-side coupling
contracts, not the staged electro-domain couplers.

Current contents:

- `electromechanicalSignalProvider.H`

The staged Purkinje, ECG, and bath-style electro couplers live under
`src/electroModels/electroCouplers/`, not here.

### `electroModels` — `libelectroModels`

The main spatial electrophysiology stack. It contains:

- top-level orchestration in `core/`, including the abstract base verifiers
  (`electroVerificationModel`, `ecgVerificationModel`, `eikonalVerificationModel`, `couplingVerificationModel`, and `graphVerificationModel`) inside `core/verificationModels/`

- domain state owners in `electroDomains/`

- numerical solver kernels in `myocardiumModels/`, `conductionSystemModels/`,

  and `ecgModels/`

- staged inter-domain couplers in `electroCouplers/`

The top-level electro entry is selected from `myocardiumSolver` in
`electroProperties`. That key first dispatches in the parent `electroModel`
runtime-selection table:

- `monodomainSolver`

- `bidomainSolver`

- `eikonalSolver`

- `singleCellSolver`

For the spatial entries, the assembled multi-domain wrapper is
`electrophysiologyModel`. It then performs the secondary dispatch into the
`myocardiumSolver` table (`monodomainSolver`, `bidomainSolver`) or builds
`eikonalMyocardiumDomain` for the canonical eikonal workflow. `singleCellSolver`
is registered directly in the parent `electroModel` table and bypasses the
myocardium-domain factory. It is not a separate `src/` library; it is compiled
inside `electroModels/myocardiumModels/`.

### `electroMechanicalModels` — `libelectroMechanicalModels`

Full electromechanical wrappers that are built only when the solids4foam
dependency is available. Lightweight EP-only builds skip this library.

## Maintained and external boundaries

- Hand-maintained project libraries live under `src/`; their source manifests
  are the adjacent `Make/files` and `Make/options` files.
- Generated ionic equation headers are outputs of the project-owned
  `applications/scripts/cellML2foam` pipeline. Change the generator, mapping,
  or template contract rather than normalizing generated equations by hand.
- `modules/physicsModel` is a project-owned compatibility layer used by the
  lightweight build.
- `modules/solids4foam` is an external submodule. cardiacFoam owns its use of
  that interface, not the submodule's implementation.

## Reading guides

### Understand the current electro runtime path

1. `src/electroModels/README.md`
1. `src/electroModels/core/README.md`
1. `src/electroModels/core/ARCHITECTURE.md`

### Trace Purkinje-to-myocardium coupling

1. `src/electroModels/electroDomains/README.md`
1. `src/electroModels/conductionSystemModels/README.md`
1. `src/electroModels/electroCouplers/README.md`

### Add or debug an ionic model

1. `src/ionicModels/README.md`
1. `src/ionicModels/IONIC_MODEL_ARCHITECTURE.md`

### Understand verification hooks

1. `src/verificationModels/README.md`
1. `src/verificationModels/VERIFICATION_MODELS_ARCHITECTURE.md`

### Get the full electrophysiology picture

1. `src/electroModels/README.md`
1. `src/electroModels/ARCHITECTURE.md`
1. `src/electroModels/core/README.md`
1. `src/electroModels/core/ARCHITECTURE.md`
1. folder-level READMEs for the subsystem you are changing
