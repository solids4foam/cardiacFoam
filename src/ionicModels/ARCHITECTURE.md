# ionicModels

`src/ionicModels` holds the runtime-selectable cellular models used by
cardiacFoam. Models are constructed through the `Foam::ionicModel` factory and
selected from the `ionicModel` entry in `electroProperties`.

The compiled model set, from `src/ionicModels/Make/files`:

Scalar models: `AlievPanfilov`, `BuenoOrovio`, `Courtemanche`, `Fabbri`,
`Gaur`, `Grandi`, `PerisYague`, `Stewart`, `TNNP`, `ToRORd_dynCl`,
`Trovato`, `TWorld`.

Batched (SoA) models: `AlievPanfilovBatched`, `BuenoOrovioBatched`,
`CourtemancheBatched`, `FabbriBatched`, `GaurBatched`, `GrandiBatched`,
`PerisYagueBatched`, `StewartBatched`, `TNNPBatched`, `ToRORd_dynClBatched`,
`TrovatoBatched`, `TWorldBatched`.

Verification models: `monodomainFDAManufactured`, `bidomainFDAManufactured`,
`bathBidomainFDAManufactured`.

## Directory layout

```text
src/ionicModels/
├── ionicModel/                    # Base class, factory, selectors, batched/GPU support headers
├── verificationModels/
│   ├── monodomainFDAManufactured/     # Manufactured monodomain ionic wrapper
│   ├── bidomainFDAManufactured/       # Manufactured bidomain ionic wrapper
│   └── bathBidomainFDAManufactured/   # Manufactured bath-bidomain ionic wrapper
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

## Base class

The central abstraction is `Foam::ionicModel` in
`src/ionicModels/ionicModel/ionicModel.H`. It provides:

- runtime selection through `ionicModel::New(...)`
- ownership of the OpenFOAM `ODESolver`
- per-integration-point ODE step sizes via `step_`
- dictionary-backed configuration storage
- tissue or dimension selection helpers
- stimulus protocol storage through `StimulusProtocol`
- generic I/O hooks for state, algebraic, rate, and constant export
- electromechanical coupling signals through `ElectromechanicalSignalProvider`

The factory's own runtime name is `TypeName("ionicModel")`. Derived classes
register with:

```cpp
addToRunTimeSelectionTable(ionicModel, MyModel, dictionary);
```

## `ionicModel/` support layers

- `ionicSelector` — shared tissue/dimension selection logic
- `ionicModelGPU` — GPU-oriented ionic-model layer
- `ionicHeterogeneity` — region-parsing and weight-math utilities (no model dependencies)
- `ionicHeterogeneityOrchestrator` — free-function orchestration layer for heterogeneity dispatch
- `configuredIonicModel` — heterogeneity support layer for scalar models (provides `HETEROGENEOUS_CONSTANTS_` and forwarding overrides)
- `batchedIonicCore` — compact SoA storage helpers for batched execution
- `batchedIonicModel` — reusable batched execution base
- `configuredBatchedIonicModel` — heterogeneity support layer for batched models (provides `HETEROGENEOUS_CONSTANTS_` and forwarding overrides)
- `batchedKernelExecution` — batched execution helpers

These headers are the shared infrastructure behind all 12 `<Model>Batched`
classes and every scalar model in `Make/files`.

## Effective derived-class contract

- `solveODE(...)`
- `derivatives(...)`
- `nEqns() const`

Common optional overrides used in the current tree:

- `supportedTissueTypes()`
- `supportedDimensions()`
- `verificationFamily()`
- `geometricDimension()`
- `ioStatesPtr()`, `ioAlgebraicPtr()`, `ioRatesPtr()`, `ioConstantsPtr()`
- `ioStateNames()`, `ioAlgebraicNames()`, `ioConstantNames()`
- `hasSignal(...)`, `signal(...)` when a model wants custom coupling behaviour

**Heterogeneity support:** scalar models gain it by inheriting from
`configuredIonicModel` rather than overriding `configureIonicHeterogeneity(...)`
directly. Exception: `ToRORd_dynCl` overrides this method to apply per-cell
initial-state blending via its own `HETEROGENEOUS_INITIAL_STATES_` member.
Batched models inherit from `configuredBatchedIonicModel` and return all three
anatomical tissue types from `supportedTissueTypes()`.

## Tissue and dimension selection

Selection is centralised in `ionicSelector`:

- `selectTissue(...)` reads `tissue`
- `selectDimension(...)` reads `dimension`

Default tissue names: `epicardialCells`, `mCells`, `endocardialCells`,
`myocyte`. Default geometric dimensions: `1D`, `2D`, `3D`. Physiological
models use `tissue`; manufactured models use `dimension`.

## ODE solver integration

The base class owns an OpenFOAM `ODESolver` and constructs it lazily:

```cpp
odeSolver_.reset(ODESolver::New(*this, dict_));
```

Ionic-model dictionaries therefore accept the standard ODE-solver keys
(`solver`, `initialODEStep`, `maxSteps`, `absTol`, `relTol`); the exact set
depends on the selected OpenFOAM ODE solver.

## Stimulus handling

The base class stores a `StimulusProtocol` loaded via
`stimulusIO::loadStimulusProtocol(dict)`. Standard ionic models typically use
`singleCellStimulus`-style inputs; manufactured models provide their own
analytical behaviour and generally skip this path. Stimulus setup lives at the
ionic-model layer, separate from the PDE-level external-stimulus blocks used
by spatial solvers.

## Generic I/O and export

`ionicModel` supports generic writing and exporting when the derived class
provides metadata and storage hooks: state names, algebraic names, constant
names, and state/algebraic/rate storage pointers. This metadata drives
single-cell trace output, export of selected variables to `volScalarField`s,
debug output filtering, and coupling-signal discovery.

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

## Coupling signals

`ionicModel` implements `ElectromechanicalSignalProvider`. The base
implementation can expose signals such as `Vm` and `Cai` from the ionic
metadata when a derived model provides recognisable state names. Some models
override `hasSignal(...)` or `signal(...)`; several simply use the base
implementation.

## Tissue heterogeneity

Ionic models can optionally support spatial heterogeneity of cellular
phenotypes through the `ionicHeterogeneity` dictionary block, letting endo/M/
epi (or other) tissue types coexist on one mesh via smooth or sharp
transitions, or named mesh regions.

- `namedRegions`: field-value-based regions with optional blending; endo/M/epi
  transmural bands are expressed as three regions tiling `[0,1]`.
- `cellZoneRegions`: hard-boundary mesh cell-zone assignment.
- `mode` is required whenever `ionicHeterogeneity` configures region-based
  heterogeneity; there is no default.
- `gradientAxes`: optional named exponential-scaling overlay axes (e.g.
  `apicobasal`), each composing multiplicatively on top of the transmural
  blend.

See `ionicHeterogeneity.H` and `ionicHeterogeneityOrchestrator.H` for the
implementation. In both scalar and batched/GPU form, only `BuenoOrovio`,
`TNNP`, `ToRORd_dynCl` and `TWorld` (and their `Batched` counterparts)
natively define three distinct tissue-type presets (`endocardialCells`,
`mCells`, `epicardialCells`), each with its own built-in constant set. The
other eight models are natively single-type — `supportedTissueTypes()`
reports `myocyte` only — but a region can still be given endo/M/epi-like
behaviour by hand: declare it with `baseline myocyte` and layer an
`ionicConstantOverrides.<name>` block on top to override the constants that
differ for that region.

## Manufactured ionic models

Three manufactured ionic-model wrappers — `monodomainFDAManufactured`,
`bidomainFDAManufactured`, `bathBidomainFDAManufactured` — differ from the
physiological models in two ways: they select by `dimension` instead of
tissue type, and they expose verification metadata through
`verificationFamily()` and `geometricDimension()`. They are used with the
verification infrastructure in `src/verificationModels`, not as standalone
physiological cell models.

## Notes on naming

Two naming conventions matter here, and should not be conflated:

- folder / runtime ionic-model name: `monodomainFDAManufactured`,
  `bidomainFDAManufactured`, `bathBidomainFDAManufactured`
- verification-model name: `manufacturedFDAMonodomainVerifier`,
  `manufacturedFDABidomainVerifier`, `manufacturedFDABathBidomainVerifier`

## Batched (SoA) support layer

Each physiological model has a `<Model>Batched` counterpart, deriving from
`configuredBatchedIonicModel` and using Structure-of-Arrays memory layout for
vectorised and GPU execution:

- `<Model>ComputeVariablesBatch` — the single `CARDIAC_HOST_DEVICE inline`
  function in `<Model>Batch.H` that processes a range `[beginCell, endCell)`.
- `<Model>RushLarsenDispatch` — a static dispatch table mapping each state
  index to its effective Rush-Larsen parameters via `rlScalarAlgAndSupport`
  or `rlNone()`.
- Compact support — a `<MODEL>_BATCH_SUPPORT_INDEX` enum defining a tight
  per-cell buffer holding effective `tau`/`gInf` pairs plus `Iion_cm`.

Runtime-selectable compact variants (`<Model>compactBatched`) are aliases that
enable compact support through `useCompactSupport_` without changing the
integration equations.
