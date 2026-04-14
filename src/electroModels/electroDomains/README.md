# electroDomains

This directory contains the domain-level state owners used by the electro
orchestration layer. Each domain implements `electroDomainInterface`, owns the
state associated with one physical region, and delegates the numerical kernel
to a runtime-selectable solver where appropriate.

## Current contents

```text
src/electroModels/electroDomains/
├── myocardiumDomain/         # Primary 3D tissue domain
├── conductionSystemDomain/   # Upstream 1D/graph conduction domain
├── ecgDomain/                # Downstream ECG evaluation domain
├── bathDomain/               # Bath-side domain code still present in tree
└── README.md
```

## Domain roles

### `MyocardiumDomain`

Defined under `myocardiumDomain/`.

- Primary domain in the `electrophysicsSystem`.
- Owns tissue fields such as `Vm`, `Iion`, `sourceField`, and
  `activationTime`.
- Implements both `tissueCouplingEndpoint` and `electroStateProvider`, so it
  can receive volumetric coupling current and expose fields to ECG or other
  consumers.
- Delegates the diffusion kernel to a runtime-selectable `myocardiumSolver`.

### `ConductionSystemDomain`

Defined under `conductionSystemDomain/`.

- Optional upstream auxiliary domain used for Purkinje or other graph-backed
  conduction networks.
- Owns the graph topology, nodal state (`Vm1D`, `Iion1D`, `activationTime`),
  PVJ metadata, and a runtime-selectable `conductionSystemSolver`.
- Implements `networkCouplingEndpoint`, exposing terminal-node voltages and
  accepting terminal coupling currents prepared by electro couplers.
- Keeps graph-specific utilities such as `conductionGraph` close to the domain
  because they are part of its state model.

### `ECGDomain`

Defined under `ecgDomain/`.

- Optional downstream domain advanced after the myocardium.
- Holds electrode configuration, ECG output, and a runtime-selectable
  `ECGSolver`.
- Consumes read-only myocardium state through `electroStateProvider`; it does
  not couple current back into the tissue.

### `bathDomain`

Bath-side domain code is still present in this folder and still compiled, but
it is not currently assembled by the active `core` orchestration path.

## Relationship to sibling directories

- Concrete tissue solvers live in
  [../myocardiumModels/README.md](../myocardiumModels/README.md).
- Concrete conduction-network solvers live in
  [../conductionSystemModels/README.md](../conductionSystemModels/README.md).
- Concrete ECG solvers live in [../ecgModels/README.md](../ecgModels/README.md).
- Inter-domain data exchange is documented in
  [../electroCouplers/README.md](../electroCouplers/README.md).
