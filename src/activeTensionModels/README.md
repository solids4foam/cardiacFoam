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
├── GoktepeKuhl/          # Goktepe-Kuhl active tension model
├── NashPanfilov/         # Nash-Panfilov active tension model
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
- Query an optional `CouplingSignalProvider` for upstream signals such as `Vm`
  and `Cai`.
- Run the base `calculateTension(...)` loop and delegate per-point work to
  `solveAtPoint(...)`.
- Provide shared I/O and export helpers through `activeTensionIO`.

## Available active-tension models

- `GoktepeKuhl`
- `NashPanfilov`

Both models select their driving electrophysiology signal from dictionary input
(`couplingSignal`, default `Vm`) and integrate with the same
`CouplingSignalProvider` interface used by `ionicModel`.
