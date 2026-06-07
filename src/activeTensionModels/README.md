# activeTensionModels library architecture

This directory provides runtime-selectable active-tension models built into
`libactiveTensionModels`.

These models sit on the electromechanics side of the stack: they consume
coupling signals from electrophysiology models and evolve active-tension state
variables for each integration point.

## Directory structure

```text
src/activeTensionModels/
├── activeTensionModel/   # Base class and runtime selection
├── GoktepeKuhl/          # Goktepe-Kuhl phenomenological active tension model
├── NashPanfilov/         # Nash-Panfilov phenomenological active tension model
├── LandNiederer/         # Land-Niederer biophysical active tension model
├── *Batched/             # GPU-ready Batched versions of the models (e.g. NashPanfilovBatched)
├── Make/
├── lnInclude/
└── README.md
```

## Core class: `Foam::activeTensionModel`

Defined in `activeTensionModel/activeTensionModel.H` and implemented in
`activeTensionModel/activeTensionModel.C`.

Main responsibilities:

- Runtime selection via `activeTensionModel::New(...)`.
- Store the active-tension dictionary and integration-point count.
- Query an optional `ElectromechanicalSignalProvider` for upstream signals such as `Vm`
  and `Cai`.
- Run the base `calculateTension(...)` loop and delegate per-point work to
  `solveAtPoint(...)`.
- Provide shared I/O and export helpers through `activeTensionIO`.

## Batched Execution (`Foam::batchedActiveTensionModel`)

For massive parallelism on CPU (OpenMP) and GPU (CUDA), models extending the `batchedActiveTensionModel` base class utilize a Structure-of-Arrays (SoA) data layout.
These batched wrappers seamlessly override the main `calculateTension(...)` loop to dispatch execution efficiently across the target backend, perfectly mirroring the `ionicModels` architecture.

## Available active-tension models

### Phenomenological

- `GoktepeKuhl` & `GoktepeKuhlBatched`
- `NashPanfilov` & `NashPanfilovBatched`

### Biophysical

- `LandNiederer` & `LandNiedererBatched`

All models select their driving electrophysiology signal from dictionary input
(`couplingSignal`, default `Vm` or `Cai`) and integrate with the
`ElectromechanicalSignalProvider` interface used by `ionicModel`.
