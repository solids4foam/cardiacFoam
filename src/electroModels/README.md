# electroModels

This folder builds `libelectroModels`, the spatial electrophysiology library.
It contains the top-level orchestration layer, domain state owners, numerical
solver kernels, and staged inter-domain couplers.

## Current structure

```text
src/electroModels/
├── core/                     # Top-level orchestration and advance schemes
├── electroDomains/           # Domain state owners
├── myocardiumModels/         # Myocardium-side solver kernels
├── conductionSystemModels/   # Purkinje/conduction solver kernels
├── ecgModels/                # ECG and bath-side solver kernels
├── electroCouplers/          # Staged electro-domain couplers
├── Make/
└── README.md

```

## Top-level runtime selection

The current top-level electro workflow is selected from:

```cpp
myocardiumSolver  monodomainSolver;

```

or:

```cpp
myocardiumSolver  bidomainSolver;

```

or:

```cpp
myocardiumSolver  eikonalSolver;

```

`electroModel::New(...)` reads that key and dispatches to the assembled
orchestration wrapper `electrophysiologyModel`.

`singleCellSolver` is also compiled in this library, but it is not part of the
multi-domain `electrophysiologyModel` path.

## Folder roles

### `core/`

Owns orchestration only:

- top-level `electroModel`

- assembled `electrophysicsSystem`

- dictionary-driven builder

- timestep advance schemes

### `electroDomains/`

Owns the long-lived state of each physical domain:

- myocardium

- conduction system / Purkinje

- ECG

- bath code is still present in the tree, but not part of the active core

  orchestration path at the moment

### `myocardiumModels/`

Contains myocardium-side solver kernels and related electro models:

- `monodomainSolver`

- `bidomainSolver`

- `eikonalSolver`

- `singleCellSolver`

### `conductionSystemModels/`

Contains Purkinje/conduction solver kernels used by
`ConductionSystemDomain`:

- `monodomain1DSolver`

- `eikonalSolver`

### `ecgModels/`

Contains downstream ECG and bath-related kernels:

- `pseudoECGSolver`

- `bathECGSolver`

- `bidomainBathECGSolver`

### `electroCouplers/`

Contains staged electro-domain coupling contracts and implementations:

- `ElectroDomainCoupler`

- endpoint interfaces

- PVJ coupling family

- `heartBathInterfaceCoupler` code remains present in the tree

## Read next

- [`ARCHITECTURE.md`](./ARCHITECTURE.md)

- [`core/README.md`](./core/README.md)

- [`core/ARCHITECTURE.md`](./core/ARCHITECTURE.md)
