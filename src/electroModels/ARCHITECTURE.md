# electroModels

Multi-domain spatial electrophysiology framework for OpenFOAM-based cardiac simulation.

Compiled as `libelectroModels`. Provides the full co-simulation infrastructure for three physically distinct cardiac subsystems — myocardium, body-surface ECG, and the conduction system — operating on separate meshes and coupled through well-defined interfaces at every timestep.

---

## Directory structure

```
electroModels/
├── core/                       Framework kernel — stable abstractions
│   ├── system/                 Domain container + dictionary-driven builder
│   ├── advanceSchemes/         Time-step orchestration strategies
│   └── electrophysiologyModel/ Concrete myocardium-centred entry point
├── electroDomains/             Physical domain implementations
│   ├── myocardiumDomain/       3D myocardium: Vm PDE + ionic ODE
│   ├── ecgDomain/              ECG body region: elliptic Vm-driven PDE
│   ├── bathDomain/             Bath-side domain (not currently wired)
│   └── conductionSystemDomain/ Purkinje/conduction network: 1D graph ODE
├── myocardiumModels/           Run-time-selectable myocardium diffusion solvers
├── ecgModels/                  Run-time-selectable ECG field solvers
├── conductionSystemModels/     Run-time-selectable conduction system solvers
└── electroCouplers/            Inter-domain coupling logic and endpoint contracts
```

---

## `core/`

The framework backbone. It owns orchestration only. Domain-family selection is pushed down into the domain layer via interfaces/factories, so `core` does not branch on myocardium or conduction solver implementations.

| File | Description |
|---|---|
| `electroModel.H/C` | OpenFOAM `physicsModel` subclass. Solver-level entry point. Owns the `electrophysicsSystem`, drives the time loop via `evolve()`, and delegates all stepping to the system. |
| `system/electrophysicsSystem.H/C` | Holds the assembled domain hierarchy: typed `myocardium_`, `conductionDomains_`, `ecgDomains_`, and coupling model lists. Implements the ordered advance → couple → write sequence. |
| `system/electrophysicsSystemBuilder.H/C` | Builds the full system from `electroProperties`. Instantiates the myocardium domain through `myocardiumDomainInterface::New(...)`, then assembles conduction and ECG-stage domains. |
| `electrophysiologyModel/` | Top-level myocardium-centered orchestration wrapper registered under `monodomainSolver`, `bidomainSolver`, and `eikonalSolver`. Owns ionic model, verification model, and field export lists. |
| `advanceSchemes/electrophysicsAdvanceScheme.H/C` | Abstract time-advance strategy. Defines when the ionic solve and the diffusion solve are called relative to each other. |
| `advanceSchemes/staggered/` | Staggered (operator-split) time advance: single-pass weak coupling. |
| `advanceSchemes/pimpleStaggered/` | PIMPLE iterative strong coupling: repeats the conduction/myocardium block until convergence. |
| `electroDomainInterface.H` | Minimal lifecycle contract that all electro domains implement. Pure virtual: `time()`, `advance(t0, dt)`. Optional no-ops: `prepareTimeStep()`, `write()`, `end()`. |
| `electroStateProvider.H` | Read-only field interface implemented by domains that expose fields upstream: `VmPtr()`, `phiEPtr()`, `conductivityPtr()`. Consumed by ECG solver and the builder. |
| `dimVoltage.H` | Shared dimension set for voltage fields (`[0 2 -3 0 0 -1 0]`, i.e. V). Avoids repeated inline dimension literals across all solvers. |
| `overrideTypeName.H` | Macro to assign a lowercase runtime type name independent of C++ class name, used with `addToRunTimeSelectionTable`. |

---

## `electroDomains/`

Each domain owns its mesh, fields, a run-time-selectable solver, and implements `electroDomainInterface`. Domains are almost entirely decoupled from each other — they communicate only through `electroStateProvider` references and `electroCouplers`.

### `myocardiumDomain/`

The primary 3D cardiac region. This folder now owns the myocardium-domain family selection.

**`myocardiumDomainInterface`**
- Inherits: `electroDomainInterface`, `tissueCouplingEndpoint`, `electroStateProvider`
- Factory: `myocardiumDomainInterface::New(...)`
- Role: selects the concrete myocardium domain from the active `myocardiumSolver` contract without leaking solver branching into `core`

**`MyocardiumDomain`**
- Inherits: `electroDomainInterface`, `tissueCouplingEndpoint`, `electroStateProvider`
- Owns: `autoPtr<myocardiumSolver>`, `ionicModel&`
- Fields: `Vm_` (transmembrane voltage), `sourceField_` (external current injection)
- On construction: detects whether the solver provides `phiE` (bidomain case) and binds it; no special-casing needed in the builder
- Domain-specific capabilities (not part of the lifecycle interface): `suggestExplicitDeltaT()`, `shouldPostProcess()`, `exportStates()`, `postProcess()`, `provider()`

**`EikonalMyocardiumDomain`**
- Inherits: `myocardiumDomainInterface`
- Owns: activation-time state (`psi`) and the 3D eikonal transport fields
- Does not require an ionic model
- Exposes activation time through the same tissue coupling endpoint used by the PVJ couplers

**`myocardiumSolver`** — abstract diffusion solver base
- `solveDiffusionExplicit(dt)` / `solveDiffusionImplicit(dt)` — spatial PDE kernel
- `phiEPtr()` — returns extracellular field pointer (or null for monodomain)
- `conductivityPtr()` — returns conductivity tensor if available

### `ecgDomain/`

The body-surface region. Solves a purely passive Laplace/Poisson equation driven by the myocardial `Vm` gradient as a source term. Has no ionic model.

**`ECGDomain`**
- Inherits: `electroDomainInterface`
- Holds: `const electroStateProvider& stateProvider_` — reference to the upstream myocardium
- Owns: `autoPtr<ecgSolver>`

**`bathDomain`**
- Bath-side domain code is still present and compiled in this tree
- It is not currently assembled by the active `core` orchestration path

**`ecgSolver`** — abstract ECG solver base. Registered implementations selected by the `ecgSolver` key in `electroProperties`.

### `conductionSystemDomain/`

The Purkinje/His-bundle network. Advances activation on a graph or 1D cable, then provides activation timing to the PVJ coupler for injection into the myocardium.

**`ConductionSystemDomain`** — concrete graph-topology conduction domain
- Inherits: `electroDomainInterface`, `networkCouplingEndpoint`
- Owns: `autoPtr<conductionSystemSolver>`, `conductionGraph`, ionic model state, PVJ metadata
- `conductionSystemSolver` — abstract 1D solver; registered implementations selected under `purkinjeNetworkModelCoeffs`

---

## `myocardiumModels/`

Concrete reaction-diffusion implementations of `myocardiumSolver`. Each registers with `addToRunTimeSelectionTable(myocardiumSolver, ...)`. The 3D eikonal myocardium path is owned by `EikonalMyocardiumDomain`, not by this solver family.

| Class | Type name | PDE / method | Notes |
|---|---|---|---|
| `monodomainSolver` | `monodomainSolver` | `∂Vm/∂t − ∇·(σᵢ∇Vm) = Iion` | Standard single-domain FVM |
| `bidomainSolver` | `bidomainSolver` | Coupled `Vm` and `phiE` | Allocates and owns `phiE` field; registered as a full factory entry |
| `eikonalSolver` | `eikonalSolver` | Activation-time wavefront `\|∇ψ\| = 1/c(x)` | No ionic ODE; anisotropy via Riemannian metric from conductivity tensor |
| `singleCellSolver` | `singleCellSolver` | ODE only: `Cm dVm/dt = −Iion + Istim` | No spatial PDE; used for ionic model validation and waveform generation |

**Monodomain PDE:**
```
Cm * ∂Vm/∂t + Iion = ∇·(σ∇Vm)/χ + Istim
```

**Bidomain PDE system:**
```
Cm * ∂Vm/∂t + Iion = ∇·(σᵢ∇Vm)/χ − ∇·(σₑ∇φₑ)/χ + Istim
0 = ∇·(σᵢ∇Vm) + ∇·(σₑ∇φₑ)   [Laplace equation for φₑ]
```

---

## `ecgModels/`

Concrete implementations of `ecgSolver`.

| Class | Type name | Method | Notes |
|---|---|---|---|
| `pseudoECGSolver` | `pseudoECG` | Volume integral of `∇Vm · r̂ / r²` | No body-conductor mesh required |
| `bathECGSolver` | `bathECGSolver` | Bath-related ECG solve | Code still present in tree |
| `bidomainBathECGSolver` | `bidomainBathECGSolver` | Bidomain + bath-side ECG solve | Code still present in tree |

---

## `conductionSystemModels/`

Concrete implementations of `conductionSystemSolver`.

| Class | Type name | Method |
|---|---|---|
| `Monodomain1DSolver` | `monodomain1DSolver` | Implicit backward-Euler cable equation + ionic ODE [default] |
| `EikonalSolver1D` | `eikonalSolver` | Eikonal fast-marching on graph — activation times only; single param `c0` [m/s] |

**Cable equation (per edge):**
```
Cm * dVm/dt + Iion = G * d²Vm/dx²  +  Istim
where G = conductance/length
```

`Monodomain1DSolver` uses the **Hines tree-elimination algorithm** for O(n) implicit solution: forward elimination leaf → root, then back-substitution root → leaf. Requires tree topology (enforced at graph load).

`EikonalSolver1D` uses a single BFS pass: `activationTime[child] = activationTime[parent] + edgeLength / c0`.

---

## `electroCouplers/`

Transfers state between domains at each timestep. Runs between domain advances in the `electrophysicsSystem` step sequence.

| Class | Role |
|---|---|
| `electroDomainCouplingEndpoints.H` | Mix-in interfaces: `tissueCouplingEndpoint` (implemented by `MyocardiumDomain`) and `networkCouplingEndpoint` (implemented by conduction system domains). Provide typed access to injection targets and activation sources. |
| `electroDomainCoupler.H/C` | Base class for all couplers. Named pair of domain references with `prepareSecondaryCoupling()`, `preparePrimaryCoupling()`, `preparePostPrimaryCoupling()` hooks. |
| `pvjCoupler/pvjMapper.H/C` | Purkinje–Ventricular Junction topology mapper. Builds the spatial map between conduction-system terminal nodes and the nearest myocardium cells. |
| `pvjCoupler/pvjCoupler.H/C` | PVJ coupling-family base. Owns the shared PVJ mapper, coupling-mode parsing, and network endpoint binding. |
| `pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.H/C` | PVJ coupling with 1D-to-3D resistance model. Reads terminal `Vm`, converts it to volumetric current, and injects it into `MyocardiumDomain::sourceField_`. |
| `pvjCoupler/eikonal/eikonalPvjCoupler.H/C` | PVJ coupling for activation-time models. Transfers Purkinje terminal activation times into the myocardium eikonal domain. |
| `heartBathInterfaceCoupler.H/C` | Bath-interface coupling code still present in the tree. |

**PVJ coupling equation** (`reactionDiffusionPvjCoupler`):
```
I_pvj = (Vm_1D − Vm_3D) / R_pvj          [A/m²]
tissue source += I_pvj  (volumetric, scattered over cells within pvjRadius)
network source -= I_pvj  (bidirectional mode only)
```

`PVJMapper` spatial algorithm: for each PVJ location, find all 3D cells within `pvjRadius`, gather tissue `Vm` as a volume-weighted average, then scatter coupling current back to those cells.

---

## Runtime selection configuration (`electroProperties`)

```cpp
myocardiumSolver  monodomainSolver;

monodomainSolverCoeffs
{
    ionicModel  BuenoOrovio;
    // ...

    conductionNetworkDomains
    {
        purkinjeNetwork
        {
            conductionSystemDomain  purkinjeNetworkModel;
            // ...
            purkinjeNetworkModelCoeffs
            {
                conductionSystemSolver  monodomain1DSolver;  // default; or: eikonalSolver
                ionicModel  BuenoOrovio;
                // ...
            }
        }
    }
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
│  2. reactionDiffusionPvjCoupler.preparePrimaryCoupling(t0, dt)│
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

## Domain field ownership

| Domain | Fields owned |
|---|---|
| `MyocardiumDomain` | `Vm_`, `Iion_`, `activationTime_`, `sourceField_`, `phiE_` (bidomain only) |
| `ConductionSystemDomain` | `Vm1D_`, `Iion1D_`, `activationTime_` (per node), `terminalCurrent_`, `terminalSource_` |
| `EikonalMyocardiumDomain` | `activationTime_` only — no ionic state |
| `ECGDomain` | electrode config, ECG output — reads Vm from myocardium via `electroStateProvider` |

---

## Inheritance summary

```
electroDomainInterface
    ├── myocardiumDomainInterface + electroStateProvider + tissueCouplingEndpoint
    │     ├── MyocardiumDomain
    │     └── EikonalMyocardiumDomain
    ├── ECGDomain
    └── ConductionSystemDomain + networkCouplingEndpoint

myocardiumSolver
    ├── monodomainSolver
    └── bidomainSolver        (owns phiE field)

electroStateProvider
    ← implemented by: MyocardiumDomain
    ← consumed by:    ECGDomain, electrophysicsSystemBuilder

ElectromechanicalSignalProvider  (from couplingModels/)
    ← implemented by: ionicModel
    ← exposed by:     electroModel::provider() → myocardium → ionicModel
```
