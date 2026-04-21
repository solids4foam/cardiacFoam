# electroModels/core

This directory is the orchestration layer for the multi-domain electrophysiology
stack. It defines the top-level model interface and coordinates the time-advance
of Purkinje, myocardium, and ECG domains within one timestep.

---

## Directory layout

```text
core/
├── electroModel.H/.C                  Top-level physicsModel façade
├── electroDomainInterface.H           Domain lifecycle contract (pure abstract)
├── electroStateProvider.H             Read-only state sharing interface
├── electroVolumeFieldDomain.H         3D FVM domain contract
├── dimVoltage.H                       Shared voltage dimension set
├── overrideTypeName.H                 Runtime type-name helper macro
│
├── system/
│   ├── electrophysicsSystem.H/.C      Domain container and advance coordinator
│   └── electrophysicsSystemBuilder.H/.C  Dictionary-driven wiring factory
│
├── advanceSchemes/
│   ├── electrophysicsAdvanceScheme.H/.C  Abstract time-advance strategy
│   ├── staggered/                        Single-pass weak coupling
│   └── pimpleStaggered/                  Iterative PIMPLE strong coupling
│
└── electrophysiologyModel/
    └── electrophysiologyModel.H/.C    Concrete myocardium-centred entry point
```

### Layout rule

> **Every self-contained logical group lives in its own subfolder.**
> **Pure abstract interfaces with no peer dependencies sit flat at `core/`.**

The flat headers (`electroDomainInterface.H`, `electroStateProvider.H`, etc.)
are included by all three subfolders and carry no dependencies on each other, so
they stay at the root level.

---

## Subdirectory responsibilities

### `system/`

Holds the two classes that assemble and coordinate the full domain graph:

- **`electrophysicsSystem`** — stores the primary myocardium domain plus
  optional conduction and ECG-stage domains and their associated couplers.
  Runs the advance → couple → write sequence by delegating to
  `electrophysicsAdvanceScheme`.
- **`electrophysicsSystemBuilder`** — constructs the domain graph and coupling
  graph from runtime dictionaries instead of hard-coded solver branches.

### `advanceSchemes/`

Runtime-selectable timestep orchestration strategies that control the ordering
of ionic solve, diffusion solve, and cross-domain coupling exchanges:

- **`staggeredElectrophysicsAdvanceScheme`** — single-pass weak coupling,
  suitable for unidirectional Purkinje-to-myocardium workflows.
- **`pimpleStaggeredElectrophysicsAdvanceScheme`** — iterative strong coupling
  using `pimpleControl`, intended for bidirectional Purkinje ↔ myocardium
  exchange.

### `electrophysiologyModel/`

Concrete `electroModel` subtype for all myocardium-centred spatial workflows.
Registers under `monodomainSolver`, `bidomainSolver`, and `eikonalSolver` in
the runtime selection table. Owns the ionic model, optional verification model,
and output field lists. Delegates spatial domain assembly to
`electrophysicsSystemBuilder`.

---

## Top-level file responsibilities

| File | Description |
|---|---|
| `electroModel.H/.C` | `physicsModel` subclass. Reads `constant/electroProperties`, owns the assembled `electrophysicsSystem`, drives the time loop via `evolve()`. |
| `electroDomainInterface.H` | Minimal lifecycle contract for all domain types: `time()`, `advance(t0,dt)`, `write()`, `end()`. |
| `electroStateProvider.H` | Read-only field interface: `VmPtr()`, `phiEPtr()`, `conductivityPtr()`. Consumed by ECG domains. |
| `electroVolumeFieldDomain.H` | Contract for 3D FVM domains: exposes `mesh()`, `VmRef()`, `Iion()`, `chi()`, `Cm()`. |
| `dimVoltage.H` | Shared `dimensionSet` for all voltage fields (V, SI). |
| `overrideTypeName.H` | Drops `override` onto OpenFOAM's `TypeName()` to suppress `-Winconsistent-missing-override`. |

---

See [ARCHITECTURE.md](./ARCHITECTURE.md) for the full timestep sequence,
interface contracts, domain order, and scheme trade-offs.
