# couplingModels architecture

This directory currently contains the coupling signal interface shared between ionic
models and downstream electromechanics code.

## Current contents

```text
src/couplingModels/
├── electromechanicsFeedbackProvider.H
├── lnInclude/
└── README.md
```

## `CouplingSignalProvider` interface

Defined in `electromechanicsFeedbackProvider.H`.

### Purpose

Provide a minimal, stable API for querying cell-level signals from electrophysiology
models without hard-coding model internals.

### Signals currently defined

- `Vm` : membrane voltage
- `Cai`: intracellular calcium

### Interface

```cpp
virtual bool hasSignal(CouplingSignal s) const = 0;
virtual scalar signal(label i, CouplingSignal s) const = 0;
```

`Foam::ionicModel` inherits from this interface.
Current base behavior:

- `Vm`: auto-discovered from model metadata or transform hook.
- `Cai`: auto-discovered from model metadata.

For active-tension models, coupling input should be chosen per model.
Example: `GoktepeKuhl` is configured with exactly one driving signal
(`couplingSignal Vm`).

## Status

This layer is intentionally small and is expected to grow as electro-mechanics
coupling paths are expanded.
