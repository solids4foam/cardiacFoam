# cardiacFoam — Source Library Architecture

This document describes the design and composition of every library in `src/`, with a detailed breakdown of `electroModels`.

---

## Top-Level Overview

```
src/
├── ionicModels/          Cell-level ODE models of cardiac electrophysiology
├── verificationModels/   Manufactured-solution and reference verifiers
├── genericWriter/        Dictionary-driven field I/O helpers
├── activeTensionModels/  Active stress/tension models for electromechanical coupling
├── couplingModels/       Electromechanical signal interfaces
├── electroModels/        Multi-domain spatial electrophysiology system (3D)
└── singleCellModel/      Standalone zero-dimensional (ODE-only) electrophysiology
```

Build order (from `src/Allwmake`):

```
genericWriter → ionicModels → verificationModels → activeTensionModels → electroModels → singleCellModel
```

Each library is independent except for the dependencies listed below.

---

## Library Descriptions

### `ionicModels` — `libionicModels`

Cell-level ionic models of cardiac action potentials. Each model is a system of ODEs governing ion channel gating variables, intracellular concentrations, and transmembrane voltage `Vm`. Models are run-time selectable via `addToRunTimeSelectionTable`.

**Available models:** AlievPanfilov, BuenoOrovio, Courtemanche, Fabbri, Gaur, Grandi, ORd, Stewart, TNNP, ToRORd_dynCl, Trovato, and GPU variants.

**Key abstractions:**
- `ionicModel` — base class; provides `computeIonicCurrent()`, `exportStates()`, ODE stepping, and the `ElectromechanicalSignalProvider` interface for active tension coupling.

**Dependencies:** none (base OpenFOAM only).

---

### `verificationModels` — `libverificationModels`

Manufactured-solution verification for monodomain and bidomain solvers. Verification is activated by an optional `verificationModel` sub-dictionary in `electroProperties` and runs pre/post-processing hooks around the ionic solve.

**Sub-libraries:**
- `electroVerification` — base class `electroVerificationModel`
- `monodomainVerification` — manufactured FDA monodomain solution
- `bidomainVerification` — manufactured FDA bidomain solution
- `ecgVerification` — ECG-specific verification tools

**Dependencies:** `ionicModels`.

---

### `genericWriter` — `libgenericWriter`

Dictionary-driven I/O helpers for writing ionic state fields, ECG signals, stimulus metadata, and post-processing fields to disk. Decouples I/O policy from solver logic — solvers declare what they produce, the writer decides when and how to write it.

**Key classes:** `ionicModelIO`, `ecgModelIO`, `stimulusIO`, `activeTensionIO`.

**Dependencies:** none.

---

### `activeTensionModels` — `libactiveTensionModels`

Active stress models used in electromechanical coupling. Compute the active tension signal `Ta` from the intracellular calcium transient supplied by the ionic model.

**Available models:** GoktepeKuhl, NashPanfilov.

**Key abstraction:** `activeTensionModel` — base class; receives a `volScalarField` of calcium and returns a `volScalarField` of active tension.

**Dependencies:** none.

---

### `couplingModels` — `libcouplingModels` *(interface only)*

Defines the `ElectromechanicalSignalProvider` abstract interface. This is the contract that ionic models implement to expose signals (`Ta`, `Vm`) to the mechanical solver. Contains no solver logic — it is a pure definition library consumed by both `ionicModels` and the solid mechanics framework.

**Dependencies:** none.

---

### `singleCellModel` — `libsingleCellModel`

A standalone zero-dimensional electrophysiology solver. Drives a single cardiac cell ionic model through a full OpenFOAM time loop without solving any spatial PDE. This is the `singleCellSolver` registered as an `electroModel` subtype.

**Purpose:**
- Isolated action potential simulations
- Ionic model validation and debugging
- Can expose `ElectromechanicalSignalProvider` for single-cell electromechanical studies

**Key class:** `SingleCellSolver` — inherits from `electroModel`; owns the time loop, `ionicModel`, stimulus, and all I/O.

**Dependencies:** `electroModels`, `ionicModels`, `verificationModels`, `genericWriter`.

---

### `electroModels` — `libelectroModels`

The multi-domain spatial electrophysiology framework. Manages the co-simulation of three spatially distinct cardiac subsystems (myocardium, ECG body, conduction system) through a shared time loop with configurable coupling.

Detailed architecture described in the next section.

**Dependencies:** `ionicModels`, `verificationModels`, `genericWriter`, `activeTensionModels`.

---

## `electroModels` — Detailed Architecture

```
src/electroModels/
├── core/
├── electroDomains/
│   ├── myocardiumDomain/
│   ├── ecgDomain/
│   └── conductionSystemDomain/
├── myocardiumModels/
├── ecgModels/
├── conductionSystemModels/
└── electroCouplers/
```

### Design principles

1. **Domain separation.** Each cardiac subsystem is an independent *domain* that implements the `electroDomainInterface` lifecycle protocol (`advance`, `write`, `end`). The system orchestrates them without knowing their internals.

2. **Solver symmetry.** Every domain owns a run-time-selectable solver. Solvers register themselves via `addToRunTimeSelectionTable` and are instantiated by the domain, never by the builder.

3. **Builder decoupling.** The `electrophysicsSystemBuilder` constructs domains by reading `electroProperties` and calling `SolverType::New(mesh, dict)`. It has no special cases per solver type.

4. **Field ownership at solver level.** Solvers that require additional fields allocate and own them (e.g. `BidomainSolver` allocates `phiE`). The domain and builder are unaware of solver-specific fields.

---

### `core/` — Framework kernel

The framework backbone. All types here are compile-time stable; nothing in `core/` depends on a specific domain or solver implementation.

| File | Role |
|---|---|
| `electroModel.H/C` | Top-level `physicsModel` subclass. Entry point for the OpenFOAM solver loop. Owns the `electrophysicsSystem` and drives the time stepping. |
| `electrophysicsSystem.H/C` | Holds the three typed domain members (`myocardium_`, `upstreamDomains_`, `downstreamDomains_`) and the coupling model list. Orchestrates the ordered advance-couple-write sequence. |
| `electrophysicsSystemBuilder.H/C` | Constructs the full `electrophysicsSystem` from `electroProperties`. Reads domain and solver types from the dictionary and delegates to domain constructors. Completely generic — no solver-specific logic. |
| `electroDomainInterface.H` | Minimal lifecycle contract for all domains. Pure virtual: `time()`, `advance()`. Optional no-ops: `prepareTimeStep()`, `write()`, `end()`. This is the element type of `PtrList<electroDomainInterface>` for the conduction system and ECG domains. |
| `electroStateProvider.H` | Read-only interface exposing field pointers (`VmPtr`, `phiEPtr`, `conductivityPtr`) from a domain to downstream consumers (e.g. ECG solver). |
| `dimVoltage.H` | Shared voltage dimensionSet constant (`[0 2 -3 0 0 -1 0]` = V). Used by all field allocations to avoid repeated dimension literals. |
| `overrideTypeName.H` | OpenFOAM macro to set a lowercase runtime type name independent of the C++ class name. |
| `schemes/` | Advance-scheme abstractions. `electrophysicsAdvanceScheme` is the base; `staggeredElectrophysicsAdvanceScheme` implements operator-splitting for explicit/implicit staggered coupling between ionic and diffusion solves. |

---

### `electroDomains/` — Domain implementations

Each domain is responsible for:
- Holding its spatial fields and geometry reference
- Owning a run-time-selectable solver instance
- Implementing `electroDomainInterface`
- Providing `electroStateProvider` data to downstream domains

#### `myocardiumDomain/`

The primary 3D myocardium region. Drives the transmembrane voltage `Vm` PDE and the ionic model ODE system in sequence.

| Class | Role |
|---|---|
| `MyocardiumDomain` | Domain class. Inherits `electroDomainInterface`, `tissueCouplingEndpoint`, `electroStateProvider`. Owns a `myocardiumSolver` and an `ionicModel`. Detects and binds `phiE` from the solver if provided (bidomain case). |
| `myocardiumSolver` | Abstract base for spatial diffusion solvers. Defines the interface for `solveDiffusionExplicit()`, `solveDiffusionImplicit()`, `phiEPtr()`, `conductivityPtr()`. Run-time selectable. |

#### `ecgDomain/`

The body-surface ECG region. Solves a purely passive elliptic PDE driven by the myocardial `Vm` as a source term. Reads from the myocardium via `electroStateProvider`.

| Class | Role |
|---|---|
| `ECGDomain` | Domain class. Inherits `electroDomainInterface`, `electroStateProvider`. Holds a reference to the upstream `electroStateProvider`. |
| `ecgSolver` | Abstract base for ECG field solvers. Run-time selectable. |

#### `conductionSystemDomain/`

The 1D Purkinje/conduction network. Advances a graph-based or 1D-mesh ODE/PDE on a network topology, then injects activation into the myocardium via coupling endpoint injection.

| Class | Role |
|---|---|
| `conductionSystemDomain` | Domain base class. Inherits `electroDomainInterface`, `networkCouplingEndpoint`. |
| `graphConductionSystemDomain` | Concrete graph-topology domain. Owns a `graphConductionSystemSolver`. |
| `graphConductionSystemSolver` | Abstract base for 1D conduction system solvers. Run-time selectable. |
| `purkinjeNetworkModel` | Graph data structure: nodes, edges, activation times. |

---

### `myocardiumModels/` — Registered myocardiumSolver implementations

All concrete diffusion solvers for the myocardium. Each registers with `addToRunTimeSelectionTable(myocardiumSolver, ...)` and is selected at runtime by the `myocardiumSolver` key in `electroProperties`.

| Solver | Type name | Description |
|---|---|---|
| `monodomainSolver` | `monodomain` | Single-domain diffusion PDE. Solves `∂Vm/∂t - ∇·(σᵢ∇Vm) = Iion`. |
| `bidomainSolver` | `bidomain` | Two-domain diffusion PDE. Solves coupled `Vm` and extracellular potential `phiE`. Allocates and owns the `phiE` field. |
| `eikonalSolver` | `eikonal` | Fast-marching activation front. Computes activation time `tau` via the eikonal equation instead of a full diffusion PDE. |

---

### `ecgModels/` — Registered ecgSolver implementations

| Solver | Type name | Description |
|---|---|---|
| `pseudoECGSolver` | `pseudoECG` | Computes a pseudo-ECG signal from the `Vm` gradient integrated over the myocardium volume. No full body-conductor PDE. |

---

### `conductionSystemModels/` — Registered graphConductionSystemSolver implementations

| Solver | Type name | Description |
|---|---|---|
| `explicit1DSolver` | `explicit1D` | Explicit ODE integration on a 1D Purkinje network graph. |
| `monodomain1DSolver` | `monodomain1D` | Monodomain diffusion on the 1D cable. |
| `eikonalSolver1D` | `eikonal1D` | Eikonal fast-marching on the 1D network. |

---

### `electroCouplers/` — Inter-domain coupling

Classes that transfer information between domains at each timestep. Couplers sit between the upstream and downstream domain advances.

| Class | Role |
|---|---|
| `electroDomainCoupler` | Base class for all couplers. Holds typed references to the two domains it bridges. |
| `electroDomainCouplingEndpoints.H` | Defines `tissueCouplingEndpoint` and `networkCouplingEndpoint` — mix-in interfaces implemented by domains that participate in coupling. |
| `pvjMapper` | Purkinje–Ventricular Junction (PVJ) topology map. Identifies which myocardium cells receive activation from which conduction-system terminal nodes. |
| `pvjResistanceCoupler` | PVJ coupling logic. Reads activation from `networkCouplingEndpoint`, applies a 1D-to-3D resistance model, and injects it as a volumetric source into the myocardium `sourceField`. |

---

### Inheritance and dependency summary

```
electroDomainInterface          ← lifecycle contract (advance/write/end)
    ├── MyocardiumDomain        also: electroStateProvider, tissueCouplingEndpoint
    ├── ECGDomain               also: electroStateProvider
    └── conductionSystemDomain  also: networkCouplingEndpoint

myocardiumSolver                ← spatial diffusion contract
    ├── monodomainSolver
    ├── bidomainSolver          (owns phiE field)
    └── eikonalSolver

electroStateProvider            ← read-only field query interface
    ← implemented by MyocardiumDomain, ECGDomain
    ← consumed by ECGDomain constructor, electrophysicsSystemBuilder

ElectromechanicalSignalProvider ← active tension / Vm signal interface
    ← implemented by ionicModel
    ← consumed by activeTensionModels, electroModel::provider()
```

---

### Runtime type selection summary

| Dictionary key | Base class | Type names |
|---|---|---|
| `myocardiumSolver` | `myocardiumSolver` | `monodomain`, `bidomain`, `eikonal` |
| `ecgSolver` | `ecgSolver` | `pseudoECG` |
| `conductionSystemSolver` | `graphConductionSystemSolver` | `explicit1D`, `monodomain1D`, `eikonal1D` |
| `electroModel` | `electroModel` | `electroActivationFoam`, `singleCellSolver` |

---

### Data flow (one timestep)

```
1. conductionSystemDomain.advance(t0, dt)
       ↓ activation pulse at terminal nodes
2. pvjResistanceCoupler.preparePrimaryCoupling(t0, dt)
       ↓ source injection into myocardium sourceField_
3. MyocardiumDomain.advance(t0, dt)
       ├── ionicModel.computeIonicCurrent()
       ├── myocardiumSolver.solveDiffusion*()
       └── phiE solve (bidomain only)
4. ECGDomain.advance(t0, dt)
       ↓ reads Vm from upstream electroStateProvider
       └── ecgSolver.solve()
5. write / end hooks on all domains
```
