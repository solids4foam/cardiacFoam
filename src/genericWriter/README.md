# genericWriter architecture

This directory contains shared helper classes used by ionic and electro models for:

- writing time-series outputs,
- mapping/exporting selected variables,
- reading and evaluating stimulus protocol dictionaries.

It is built as `libgenericWriter`.

## Directory structure

```text
src/genericWriter/
├── ionicModelIO.H/.C
├── ionicVariableCompatibility.H/.C
├── stimulusIO.H/.C
├── Make/
└── README.md
```

## Components

### `ionicModelIO`

Central utility for ionic model outputs:

- header generation (`writeHeader`, `writeSelectedHeader`),
- data row writes (`write`, `writeSelected`),
- timestep write gating (`shouldWriteStep`),
- selected-name filtering and caching,
- state/algebraic/rate export to `volScalarField`,
- debug-print helpers,
- sweep output helpers (`writeSweepHeader`, `writeOneSweepRow`).

It also provides small cache structs to avoid repeating expensive variable resolution.

### `ionicVariableCompatibility`

Name compatibility and resolution layer:

- relaxed matching for Vm aliases,
- centralized `Vm`/`Cai` state-index discovery used by coupling signals,
- centralized activation-like state discovery (`Act`) used by coupling signals,
- relaxed matching for `RATES_<state>` names,
- mapping user-requested names to state/algebraic/rate indices.

Used by `ionicModelIO` when filtering exported/debug variable lists.

### `stimulusIO`

Dictionary parser and stimulus evaluator:

- `StimulusProtocol` (single-cell S1/S2 timing protocol),
- `MonodomainStimulusProtocol` (spatial box-based monodomain stimulus list),
- protocol load functions from OpenFOAM dictionaries,
- pulse evaluation (`computeStimulus`),
- helper checks (`hasActiveStimulus`, `protocolSuffix`).

## Build target

`Make/files` builds:

- `ionicModelIO.C`
- `ionicVariableCompatibility.C`
- `stimulusIO.C`

into:

- `$(FOAM_USER_LIBBIN)/libgenericWriter`
