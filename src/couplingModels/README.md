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
- `Act`: activation signal
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
- `Act` (non-manufactured models): resolved from existing states using
  priority `Act`-like state -> `Cai` -> `Vm`.

Reduced models can still provide explicit `Act` overrides.

For active-tension models, coupling input should be chosen per model.
Example: `GoktepeKuhl` is configured with exactly one driving signal
(`couplingSignal Act` or `couplingSignal Vm`), not both simultaneously.

## Status

This layer is intentionally small and is expected to grow as electro-mechanics
coupling paths are expanded.
