# electroDomains

This directory contains the domain-level state owners used by the electro
orchestration layer. Each domain implements `electroDomainInterface`, owns the
state associated with one physical region, and delegates the numerical kernel
to a runtime-selectable solver where appropriate.

## Current contents

```text
src/electroModels/electroDomains/
├── myocardiumDomain/                # Primary 3D tissue domain
├── conductionSystemDomain/          # Upstream 1D/graph conduction domain
├── ecgDomain/                       # Downstream ECG evaluation domain
├── extracellularPotentialDomain/    # Unified global phiE (heart+bath)
└── README.md
```

## Domain roles

### `myocardiumDomain`

Defined under `myocardiumDomain/`.

- Primary domain in the `electrophysicsSystem`.
- Owns tissue fields such as `Vm`, `Iion`, `sourceField`, and
  `activationTime`.
- Implements both `tissueCouplingEndpoint` and `electroStateProvider`, so it
  can receive volumetric coupling current and expose fields to ECG or other
  consumers.
- Delegates the diffusion kernel to a runtime-selectable `myocardiumSolver`.

### `conductionSystemDomain`

Defined under `conductionSystemDomain/`.

- Optional upstream auxiliary domain used for Purkinje or other graph-backed
  conduction networks.
- Owns the graph topology, nodal state (`Vm1D`, `Iion1D`, `activationTime`),
  PVJ metadata, and a runtime-selectable `conductionSystemSolver`.
- Implements `networkCouplingEndpoint`, exposing terminal-node voltages and
  accepting terminal coupling currents prepared by electro couplers.
- Keeps graph-specific utilities such as `conductionGraph` close to the domain
  because they are part of its state model.

### `ecgDomain`

Defined under `ecgDomain/`.

- Optional downstream domain advanced after the myocardium.
- Holds electrode configuration, ECG output, and a runtime-selectable
  `ecgSolver`.
- Consumes read-only myocardium state through `electroStateProvider`; it does
  not couple current back into the tissue.

### `extracellularPotentialDomain`

Defined under `extracellularPotentialDomain/`.

- Implements `electroStateDomain` — both advances in time and exposes
  read-only state (phiE, conductivities).
- Owns the global extracellular potential `phiE` solved on the union mesh of
  heart + bath cell zones via an elliptic FVM laplacian.
- Scatters heart `Vm` from the bidomain solver to the base mesh and binds a
  restricted local view of `phiE` back into the myocardium via
  `bindExternalPhiE`. With this binding the bidomain solver no longer solves
  its own local phiE.
- Provides the `electroStateProvider` accessed by `bathECGProbe`-class ECG
  domains for electrode sampling on the full union mesh.

## Relationship to sibling directories

- Concrete tissue solvers live in
  [../myocardiumModels/README.md](../myocardiumModels/README.md).
- Concrete conduction-network solvers live in
  [../conductionSystemModels/README.md](../conductionSystemModels/README.md).
- Concrete ECG solvers live in [../ecgModels/README.md](../ecgModels/README.md).
- Inter-domain data exchange is documented in
  [../electroCouplers/README.md](../electroCouplers/README.md).
