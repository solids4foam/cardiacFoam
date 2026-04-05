# electroModels

Multi-domain spatial electrophysiology framework for OpenFOAM-based cardiac simulation.

Compiled as `libelectroModels`. Provides the full co-simulation infrastructure for three physically distinct cardiac subsystems — myocardium, body-surface ECG, and the conduction system — operating on separate meshes and coupled through well-defined interfaces at every timestep.

---

## Directory structure

```
electroModels/
├── core/                       Framework kernel — stable abstractions
├── electroDomains/             Physical domain implementations
│   ├── myocardiumDomain/       3D myocardium: Vm PDE + ionic ODE
│   ├── ecgDomain/              ECG body region: elliptic Vm-driven PDE
│   └── conductionSystemDomain/ Purkinje/conduction network: 1D graph ODE
├── myocardiumModels/           Run-time-selectable myocardium diffusion solvers
├── ecgModels/                  Run-time-selectable ECG field solvers
├── conductionSystemModels/     Run-time-selectable conduction system solvers
└── electroCouplers/            Inter-domain coupling logic and endpoint contracts
```

---

## `core/`

The framework backbone. Nothing here depends on any specific domain or solver.

| File | Description |
|---|---|
| `electroModel.H/C` | OpenFOAM `physicsModel` subclass. Solver-level entry point. Owns the `electrophysicsSystem`, drives the time loop via `evolve()`, and delegates all stepping to the system. |
| `electrophysicsSystem.H/C` | Holds the assembled domain hierarchy: typed `myocardium_` (`autoPtr<MyocardiumDomain>`), `upstreamDomains_` and `downstreamDomains_` (`PtrList<electroDomainInterface>`), and coupling model lists. Implements the ordered advance → couple → write sequence. |
| `electrophysicsSystemBuilder.H/C` | Builds the full system from `electroProperties`. Reads domain and solver type names; instantiates domains using their constructors and `SolverType::New()`. Contains no solver-specific branches. |
| `electroDomainInterface.H` | Minimal lifecycle contract that all electro domains implement. Pure virtual: `time()`, `advance(t0, dt)`. Optional no-ops: `prepareTimeStep()`, `write()`, `end()`. This is the base type stored in `PtrList<electroDomainInterface>`. |
| `electroStateProvider.H` | Read-only field interface implemented by domains that expose fields upstream: `VmPtr()`, `phiEPtr()`, `conductivityPtr()`. Consumed by ECG solver and the builder. |
| `dimVoltage.H` | Shared dimension set for voltage fields (`[0 2 -3 0 0 -1 0]`, i.e. V). Avoids repeated inline dimension literals across all solvers. |
| `overrideTypeName.H` | Macro to assign a lowercase runtime type name independent of C++ class name, used with `addToRunTimeSelectionTable`. |
| `schemes/electrophysicsAdvanceScheme.H/C` | Abstract time-advance strategy. Defines when the ionic solve and the diffusion solve are called relative to each other. |
| `schemes/staggered/staggeredElectrophysicsAdvanceScheme.H/C` | Staggered (operator-split) time advance: ionic half-step, diffusion step, ionic half-step. Controls nodal stability for explicit–implicit coupling. |

---

## `electroDomains/`

Each domain owns its mesh, fields, a run-time-selectable solver, and implements `electroDomainInterface`. Domains are almost entirely decoupled from each other — they communicate only through `electroStateProvider` references and `electroCouplers`.

### `myocardiumDomain/`

The primary 3D cardiac region. Solves the transmembrane voltage PDE coupled to the ionic model ODE system at each cell.

**`MyocardiumDomain`**
- Inherits: `electroDomainInterface`, `tissueCouplingEndpoint`, `electroStateProvider`
- Owns: `autoPtr<myocardiumSolver>`, `ionicModel&`
- Fields: `Vm_` (transmembrane voltage), `sourceField_` (external current injection)
- On construction: detects whether the solver provides `phiE` (bidomain case) and binds it; no special-casing needed in the builder
- Domain-specific capabilities (not part of the lifecycle interface): `suggestExplicitDeltaT()`, `shouldPostProcess()`, `exportStates()`, `postProcess()`, `provider()`

**`myocardiumSolver`** — abstract diffusion solver base
- `solveDiffusionExplicit(dt)` / `solveDiffusionImplicit(dt)` — spatial PDE kernel
- `phiEPtr()` — returns extracellular field pointer (or null for monodomain)
- `conductivityPtr()` — returns conductivity tensor if available

### `ecgDomain/`

The body-surface region. Solves a purely passive Laplace/Poisson equation driven by the myocardial `Vm` gradient as a source term. Has no ionic model.

**`ECGDomain`**
- Inherits: `electroDomainInterface`, `electroStateProvider`
- Holds: `const electroStateProvider& stateProvider_` — reference to the upstream myocardium
- Owns: `autoPtr<ecgSolver>`

**`ecgSolver`** — abstract ECG solver base. Registered implementations selected by the `ecgSolver` key in `electroProperties`.

### `conductionSystemDomain/`

The Purkinje/His-bundle network. Advances activation on a graph or 1D cable, then provides activation timing to the PVJ coupler for injection into the myocardium.

**`conductionSystemDomain`** — abstract domain base
- Inherits: `electroDomainInterface`, `networkCouplingEndpoint`

**`graphConductionSystemDomain`** — concrete graph-topology domain
- Owns: `autoPtr<graphConductionSystemSolver>`, `purkinjeNetworkModel`
- `graphConductionSystemSolver` — abstract solver; registered implementations selected by key

**`purkinjeNetworkModel`** — graph data structure: nodes, conduction edges, activation time map.

---

## `myocardiumModels/`

Concrete implementations of `myocardiumSolver`. Each registers with `addToRunTimeSelectionTable(myocardiumSolver, ...)`.

| Class | Type name | PDE / method | Notes |
|---|---|---|---|
| `monodomainSolver` | `monodomain` | `∂Vm/∂t − ∇·(σᵢ∇Vm) = Iion` | Standard single-domain FVM |
| `bidomainSolver` | `bidomain` | Coupled `Vm` and `phiE` | Allocates and owns `phiE` field; registered as a full factory entry |
| `eikonalSolver` | `eikonal` | Fast-marching activation time `τ` | No time-integration of ion channels; for propagation studies |

---

## `ecgModels/`

Concrete implementations of `ecgSolver`.

| Class | Type name | Method | Notes |
|---|---|---|---|
| `pseudoECGSolver` | `pseudoECG` | Volume integral of `∇Vm · r̂ / r²` | No body-conductor mesh required |

---

## `conductionSystemModels/`

Concrete implementations of `graphConductionSystemSolver`.

| Class | Type name | Method |
|---|---|---|
| `explicit1DSolver` | `explicit1D` | Explicit Euler ODE on cable segments |
| `monodomain1DSolver` | `monodomain1D` | Monodomain diffusion on 1D cable |
| `eikonalSolver1D` | `eikonal1D` | Eikonal fast-marching on graph topology |

---

## `electroCouplers/`

Transfers state between domains at each timestep. Runs between domain advances in the `electrophysicsSystem` step sequence.

| Class | Role |
|---|---|
| `electroDomainCouplingEndpoints.H` | Mix-in interfaces: `tissueCouplingEndpoint` (implemented by `MyocardiumDomain`) and `networkCouplingEndpoint` (implemented by conduction system domains). Provide typed access to injection targets and activation sources. |
| `electroDomainCoupler.H/C` | Base class for all couplers. Named pair of domain references with `prepareSecondaryCoupling()`, `preparePrimaryCoupling()`, `preparePostPrimaryCoupling()` hooks. |
| `pvjMapper.H/C` | Purkinje–Ventricular Junction topology mapper. Builds the spatial map between conduction-system terminal nodes and the nearest myocardium cells. |
| `pvjResistanceCoupler.H/C` | PVJ coupling with 1D-to-3D resistance model. Reads activation times from `networkCouplingEndpoint`, converts to volumetric current, injects into `MyocardiumDomain::sourceField_`. |

---

## Runtime selection configuration (`electroProperties`)

```cpp
electroModel        electroActivationFoam;

myocardiumSolverCoeffs
{
    myocardiumSolver    monodomain;   // or: bidomain, eikonal
    ionicModel          BuenoOrovio;
}

ecgSolverCoeffs
{
    ecgSolver           pseudoECG;
}

conductionSystemSolverCoeffs
{
    conductionSystemSolver  explicit1D;
}
```

---

## Timestep data flow

```
┌───────────────────────────────────────────────────────────────┐
│                    electrophysicsSystem                        │
│                                                               │
│  1. conductionSystemDomain.advance(t0, dt)                    │
│         → activates terminal nodes of Purkinje graph          │
│                                                               │
│  2. pvjResistanceCoupler.preparePrimaryCoupling(t0, dt)       │
│         → injects current into myocardium sourceField_        │
│                                                               │
│  3. MyocardiumDomain.advance(t0, dt)                          │
│       ├─ ionicModel.computeIonicCurrent()  [per cell, ODE]    │
│       ├─ myocardiumSolver.solveDiffusion*() [Vm PDE, FVM]     │
│       └─ phiE solve (bidomain only)                           │
│                                                               │
│  4. ECGDomain.advance(t0, dt)                                 │
│       └─ ecgSolver.solve()  [reads Vm from myocardium]        │
│                                                               │
│  5. write() / end() on all domains and couplers               │
└───────────────────────────────────────────────────────────────┘
```

---

## Inheritance summary

```
electroDomainInterface
    ├── MyocardiumDomain      + electroStateProvider + tissueCouplingEndpoint
    ├── ECGDomain             + electroStateProvider
    └── conductionSystemDomain + networkCouplingEndpoint
          └── graphConductionSystemDomain

myocardiumSolver
    ├── monodomainSolver
    ├── bidomainSolver        (owns phiE field)
    └── eikonalSolver

electroStateProvider
    ← implemented by: MyocardiumDomain, ECGDomain
    ← consumed by:    ECGDomain, electrophysicsSystemBuilder

ElectromechanicalSignalProvider  (from couplingModels/)
    ← implemented by: ionicModel
    ← exposed by:     electroModel::provider() → myocardium → ionicModel
```
