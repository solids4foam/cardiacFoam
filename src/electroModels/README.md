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

The top-level electro workflow is selected via `myocardiumSolver` in `constant/electroProperties`:

```cpp
myocardiumSolver  monodomainSolver;  // or: bidomainSolver | eikonalSolver
```

`electroModel::New(...)` reads that key and dispatches to `electrophysiologyModel`.

`singleCellSolver` is compiled in this library but is not part of the multi-domain `electrophysiologyModel` path.

## Selection and ownership

```text
physicsModel -> electroModel -> electrophysiologyModel
                                  -> myocardiumDomain -> myocardiumSolver
                                  -> conduction/ECG/bath domains and couplers
```

`constant/physicsProperties` selects the top-level `electroModel`.
`constant/electroProperties` then supplies the canonical `myocardiumSolver`
runtime name and its matching `<name>Coeffs` dictionary. The spatial names
`monodomainSolver`, `bidomainSolver`, and `eikonalSolver` select the common
`electrophysiologyModel`; the builder then creates the corresponding myocardium
domain and optional `conductionNetworkDomains`, `ecgDomains`,
`bathPotentialDomain`, and `domainCouplings` entries.

The orchestration layer controls ordering but does not own numerical fields.
Domains own long-lived state and meshes; solver classes implement domain-local
numerical kernels; couplers transfer state through typed endpoints. These EP
layers build in both full and lightweight modes. Electromechanical wrappers are
in `src/electroMechanicalModels` and require full solids4foam mode.

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

- extracellular potential / bath ECG through `extracellularPotentialDomain`
  and `torsoECG`

### `myocardiumModels/`

Contains myocardium-side solver kernels and related electro models:

- `monodomainSolver`

- `bidomainSolver`

- `eikonalSolver`

- `singleCellSolver`

### `conductionSystemModels/`

Contains Purkinje/conduction solver kernels used by
`conductionSystemDomain`:

- `monodomain1DSolver`

- `eikonalSolver`

### `ecgModels/`

Contains downstream ECG kernels:

- `pseudoECG` (implemented by class `pseudoECGSolver`)

- `torsoECG` — electrode sampler on the global phiE from
  `extracellularPotentialDomain`

### `electroCouplers/`

Contains staged electro-domain coupling contracts and implementations:

- `electroDomainCoupler`

- endpoint interfaces

- PVJ coupling family

Bath coupling is handled directly by `extracellularPotentialDomain`, which owns
the global `phiE` solve and binds a restricted phiE view into the bidomain
myocardium solver. No dedicated coupler class is needed.

## Read next

- [`ARCHITECTURE.md`](./ARCHITECTURE.md)

- [`core/README.md`](./core/README.md)

- [`core/ARCHITECTURE.md`](./core/ARCHITECTURE.md)
