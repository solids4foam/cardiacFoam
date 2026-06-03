# Ionic Models Architecture

## Overview

`src/ionicModels` contains the runtime-selectable cellular models used by
`cardiacFoam`. These models are constructed through the `Foam::ionicModel`
factory and are selected from the `ionicModel` entry in `electroProperties`.

The current compiled model set is defined by
`src/ionicModels/Make/files`. In this repository that list is:

Scalar models: `AlievPanfilov`, `BuenoOrovio`, `Courtemanche`, `Fabbri`,
`Gaur`, `Grandi`, `PerisYague`, `Stewart`, `TNNP`, `ToRORd_dynCl`,
`Trovato`, `TWorld`.

Batched (SoA) models: `AlievPanfilovBatched`, `BuenoOrovioBatched`,
`CourtemancheBatched`, `FabbriBatched`, `GaurBatched`, `GrandiBatched`,
`PerisYagueBatched`, `StewartBatched`, `TNNPBatched`, `ToRORd_dynClBatched`,
`TrovatoBatched`.

Verification models: `monodomainFDAManufactured`, `bidomainFDAManufactured`,
`bathBidomainFDAManufactured`.

This file describes the architecture that actually exists in this tree. It is
not intended as a generic survey of cardiac ionic models.

## Directory Layout

```text
src/ionicModels/
├── ionicModel/                    # Base class, factory, selectors, batched/GPU support headers
├── monodomainFDAManufactured/     # Manufactured monodomain ionic wrapper
├── bidomainFDAManufactured/       # Manufactured bidomain ionic wrapper
├── bathBidomainFDAManufactured/   # Manufactured bath-bidomain ionic wrapper
├── AlievPanfilov/
├── AlievPanfilovBatched/
├── BuenoOrovio/
├── BuenoOrovioBatched/
├── Courtemanche/
├── CourtemancheBatched/
├── Fabbri/
├── FabbriBatched/
├── Gaur/
├── GaurBatched/
├── Grandi/
├── GrandiBatched/
├── PerisYague/
├── PerisYagueBatched/
├── Stewart/
├── StewartBatched/
├── TNNP/
├── TNNPBatched/
├── ToRORd_dynCl/
├── ToRORd_dynClBatched/
├── Trovato/
├── TrovatoBatched/
├── TWorld/
├── Make/
└── lnInclude/
```

## Base Class

The central abstraction is `Foam::ionicModel` in
`src/ionicModels/ionicModel/ionicModel.H`.

The base class provides:

- runtime selection through `ionicModel::New(...)`
- ownership of the OpenFOAM `ODESolver`
- per-integration-point ODE step sizes via `step_`
- dictionary-backed configuration storage
- tissue or dimension selection helpers
- stimulus protocol storage through `StimulusProtocol`
- generic I/O hooks for state, algebraic, rate, and constant export
- electromechanical coupling signals through `ElectromechanicalSignalProvider`

The runtime name for the factory itself is:

```cpp
TypeName("ionicModel");
```

Derived classes register with:

```cpp
addToRunTimeSelectionTable(ionicModel, MyModel, dictionary);
```

## `ionicModel/` support layers

The `ionicModel/` folder now contains more than the classic base class and
factory code. It also includes:

- `ionicSelector`
  shared tissue/dimension selection logic
- `ionicModelGPU`
  GPU-oriented ionic-model layer
- `batchedIonicCore`
  compact SoA storage helpers for batched execution
- `batchedIonicModel`
  reusable batched execution base
- `configuredBatchedIonicModel`
  shared configuration/metadata layer for batched models
- `batchedKernelExecution`
  batched execution helpers

These headers are support infrastructure shared by all 11 `<Model>Batched`
classes that are compiled via `Make/files`.

## Effective Derived-Class Contract

This codebase does not use the old `advance()/nGates()/initialVm()` style
interface. The effective contract in this repository is the one implemented by
the current derived classes:

- `solveODE(...)`
- `derivatives(...)`
- `nEqns() const`

Common optional overrides used in the current tree are:

- `supportedTissueTypes()`
- `supportedDimensions()`
- `verificationFamily()`
- `geometricDimension()`
- `ioStatesPtr()`, `ioAlgebraicPtr()`, `ioRatesPtr()`, `ioConstantsPtr()`
- `ioStateNames()`, `ioAlgebraicNames()`, `ioConstantNames()`
- `configureIonicHeterogeneity(const scalarField&, const dictionary&)` — applies transmural spatial heterogeneity to ionic properties (see section below)
- `hasSignal(...)`, `signal(...)` when a model wants custom coupling behavior

## Tissue And Dimension Selection

Selection is centralized in `ionicSelector`:

- `selectTissue(...)` reads `tissue`
- `selectDimension(...)` reads `dimension`

Default tissue names accepted by the base class are:

- `epicardialCells`
- `mCells`
- `endocardialCells`
- `myocyte`

Default geometric dimension names are:

- `1D`
- `2D`
- `3D`

Most physiological ionic models use `tissue`. The manufactured models use
`dimension`.

## ODE Solver Integration

The base class owns an OpenFOAM `ODESolver` and constructs it lazily from the
model dictionary:

```cpp
odeSolver_.reset(ODESolver::New(*this, dict_));
```

That means ionic-model dictionaries can pass through standard ODE-solver keys
such as:

- `solver`
- `initialODEStep`
- `maxSteps`
- `absTol`
- `relTol`

The exact accepted keys depend on the selected OpenFOAM ODE solver.

## Stimulus Handling

The base class stores a `StimulusProtocol` loaded from the model dictionary via
`stimulusIO::loadStimulusProtocol(dict)`.

In practice:

- standard ionic models typically use `singleCellStimulus`-style inputs
- manufactured models provide their own analytical behavior and generally do
  not rely on the same pacing path

Stimulus setup lives at the ionic-model layer, separate from PDE-level external
stimulus blocks used by spatial solvers.

## Generic I/O And Export

`ionicModel` supports generic writing and exporting when the derived class
provides metadata and storage hooks. The important hooks are:

- state names
- algebraic names
- constant names
- state/algebraic/rate storage pointers

This metadata is used for:

- writing single-cell traces
- exporting selected ionic variables to `volScalarField`s
- debug output filtering
- coupling-signal discovery

The configuration shape consumed by the base class is:

```text
outputVariables
{
    ionic
    {
        export (...);
        debug  (...);
    }
}
```

## Coupling Signals

`ionicModel` also implements `ElectromechanicalSignalProvider`.

The base implementation can expose signals such as `Vm` and `Cai` from the
ionic metadata when a derived model provides recognizable state names. Some
models override `hasSignal(...)` or `signal(...)`, but several simply delegate
to the base implementation.

## Tissue Heterogeneity

Ionic models can optionally support spatial transmural heterogeneity of cellular
phenotypes through the `configureIonicHeterogeneity(...)` virtual method. This
allows different tissue types (endocardial, mid-myocardial, epicardial) to coexist
within a single spatial mesh via smooth or sharp transitions.

**Configuration:** The `ionicHeterogeneity` dictionary block is read from the
model coefficients and specifies a transmural distance field, transition boundaries,
and smoothing parameters. See the detailed schema in `src/ionicModels/README.md`.

**Current support:**

- **Scalar CPU:** `BuenoOrovio` (transmural bands with smooth transitions)
- **Batched/GPU:** `BuenoOrovioBatched`, `TNNPBatched`, `TWorldBatched`, `ToRORd_dynClBatched`

## Manufactured Ionic Models

This repository currently has two manufactured ionic-model wrappers:

- `monodomainFDAManufactured`
- `bidomainFDAManufactured`

They differ from the standard physiological models in two ways:

- they select by `dimension` instead of tissue type
- they expose verification metadata through `verificationFamily()` and
  `geometricDimension()`

These models are used together with the verification infrastructure in
`src/verificationModels`, not as standalone physiological cell models.

## Notes On Naming

Two naming conventions matter here:

- folder / runtime ionic-model name:
  `monodomainFDAManufactured`, `bidomainFDAManufactured`
- verification-model name:
  `manufacturedFDAMonodomainVerifier`,
  `manufacturedFDABidomainVerifier`

These are different layers and should not be conflated.

## Batched (SoA) Support Layer

Each physiological model has a `<Model>Batched` counterpart. These derive from
`configuredBatchedIonicModel` and use Structure-of-Arrays (SoA) memory layout
for efficient vectorised and GPU execution.

The batched path uses:

- `<Model>ComputeVariablesBatch` — the single `CARDIAC_HOST_DEVICE inline` function
  in `<Model>Batch.H` that processes a range `[beginCell, endCell)`.
- `<Model>RushLarsenDispatch` — a static dispatch table mapping each state index
  to its effective RL parameters via `rlScalarAlgAndSupport` or `rlNone()`.
- Compact support — a `<MODEL>_BATCH_SUPPORT_INDEX` enum defining a tight
  per-cell buffer holding effective `tau`/`gInf` pairs plus `Iion_cm`.

Runtime-selectable compact variants (`<Model>compactBatched`) are aliases that
enable compact support through `useCompactSupport_` without changing the
integration equations.

## Deliberate Non-Claims

This document intentionally does not claim the following, because they are not
the current truth of this tree:

- a `Mitchell` ionic model
- a `TenTusscher` folder/runtime type separate from `TNNP`
- an `ORd` folder or runtime type separate from `ToRORd_dynCl`
- a base-class contract centered on `advance()` / `initialGates()`

If those features are added later, this document should be updated from the
actual source and `Make/files`, not from external expectations.
