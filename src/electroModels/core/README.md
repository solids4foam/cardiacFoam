# electroModels core

This directory contains the orchestration layer for the multi-domain
electrophysiology stack. It defines the top-level model interface, the
assembled domain system, and the time-advance strategies that coordinate
Purkinje, myocardium, and ECG domains within one timestep.

## Current contents

```text
src/electroModels/core/
├── electroModel.{H,C}                 # Top-level physicsModel facade
├── electrophysicsSystem.{H,C}         # Assembled domain container
├── electrophysicsSystemBuilder.{H,C}  # Dictionary-driven wiring
├── electroDomainInterface.H           # Shared domain lifecycle contract
├── electroStateProvider.H             # Read-only state sharing interface
├── dimVoltage.H                       # Shared voltage dimensions
├── overrideTypeName.H                 # Runtime type-name helper macro
├── schemes/
│   ├── electrophysicsAdvanceScheme.{H,C}
│   ├── staggered/                     # Single-pass weak coupling
│   └── pimpleStaggered/               # Iterative strong coupling
├── electroActivationFoam/             # Legacy electroModel-style driver
├── ELECTROMODEL_ORCHESTRATION.md      # Detailed design note
└── README.md
```

## Main responsibilities

- `electroModel` reads `constant/electroProperties`, owns the assembled
  `electrophysicsSystem`, and exposes the solver-facing `evolve()` entry point.
- `electrophysicsSystem` stores the primary myocardium domain plus optional
  upstream and downstream domains and their associated couplers.
- `electrophysicsSystemBuilder` constructs the domain graph and coupling graph
  from runtime dictionaries instead of hard-coded solver branches.
- `electroDomainInterface` defines the common domain lifecycle:
  `prepareTimeStep()`, `advance()`, `write()`, and `end()`.
- `electroStateProvider` exposes read-only fields such as `Vm`, `phiE`, and
  conductivity so one domain can drive another without tight coupling.

## Advance schemes

The `schemes/` subdirectory contains runtime-selectable timestep orchestration:

- `staggeredElectrophysicsAdvanceScheme`: single-pass weak coupling, suitable
  for unidirectional Purkinje-to-myocardium workflows.
- `pimpleStaggeredElectrophysicsAdvanceScheme`: iterative strong coupling using
  `pimpleControl`, intended for bidirectional Purkinje-myocardium exchange.

See [ELECTROMODEL_ORCHESTRATION.md](./ELECTROMODEL_ORCHESTRATION.md) for the
full timestep sequence, domain order, and scheme tradeoffs.
