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
│   ├── myocardiumDomain/                3D myocardium: Vm PDE + ionic ODE
│   ├── ecgDomain/                       ECG body region: elliptic Vm-driven PDE
│   ├── extracellularPotentialDomain/    Unified global phiE (heart+bath)
│   └── conductionSystemDomain/          Purkinje/conduction network: 1D graph ODE
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
| `dimVoltage.H` | Shared dimension set for voltage fields (`[1 2 -3 0 0 -1 0]`, i.e. V). Avoids repeated inline dimension literals across all solvers. |
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

**`myocardiumDomain`**

- Inherits: `electroDomainInterface`, `tissueCouplingEndpoint`, `electroStateProvider`

- Owns: `autoPtr<myocardiumSolver>`, `ionicModel&`

- Fields: `Vm_` (transmembrane voltage), `sourceField_` (external current injection)

- On construction: detects whether the solver provides `phiE` (bidomain case) and binds it; no special-casing needed in the builder

- Domain-specific capabilities (not part of the lifecycle interface): `suggestExplicitDeltaT()`, `shouldPostProcess()`, `exportStates()`, `postProcess()`, `provider()`

**`eikonalMyocardiumDomain`**

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

**`ecgDomain`**

- Inherits: `electroDomainInterface`

- Holds: `const electroStateProvider& stateProvider_` — reference to the upstream myocardium

- Owns: `autoPtr<ecgSolver>`

**`extracellularPotentialDomain`** — unified global extracellular potential
domain. Implements `electroStateDomain`.

- Inherits: `electroDomainInterface`, `electroStateProvider`

- Owns: `phiE` field on a base mesh covering the heart cell zone plus one or
  more bath cell zones.

- Solves an elliptic Poisson equation for `phiE` driven by the scattered heart
  `Vm`, with Dirichlet `groundPatches` and Neumann `surfaceCurrentPatches`.

- Binds a restricted local view of `phiE` back into the bidomain myocardium
  solver via `bindExternalPhiE`, so the bidomain solver no longer solves a
  local phiE.

- Exposes `phiE` to `torsoECG`-class ECG solvers for electrode sampling.

- Writes `phiE`, `sigmaTotal`, and `VmGlobal` through the normal electro-model
  output path at OpenFOAM output times.

**`ecgSolver`** — abstract ECG solver base. Registered implementations selected by the `ecgSolver` key in `electroProperties`.

### `conductionSystemDomain/`

The Purkinje/His-bundle network. Advances activation on a graph or 1D cable, then provides activation timing to the PVJ coupler for injection into the myocardium.

**`conductionSystemDomain`** — concrete graph-topology conduction domain

- Inherits: `electroDomainInterface`, `networkCouplingEndpoint`

- Owns: `autoPtr<conductionSystemSolver>`, `conductionGraph`, ionic model state, PVJ metadata

- `conductionSystemSolver` — abstract 1D solver; registered implementations selected under `purkinjeGraphModelCoeffs`

---

## `myocardiumModels/`

Concrete reaction-diffusion implementations of `myocardiumSolver`. Only `monodomainSolver` and `bidomainSolver` register in this table. The `myocardiumSolver` dictionary key is a top-level selector: it first dispatches through the parent `electroModel` table, then `electrophysiologyModel` builds the concrete myocardium domain through `myocardiumDomainInterface::New(...)`.

| Class | Type name | PDE / method | Notes |
|---|---|---|---|
| `monodomainSolver` | `monodomainSolver` | `∂Vm/∂t − ∇·(σᵢ∇Vm) = Iion` | Standard single-domain FVM |
| `bidomainSolver` | `bidomainSolver` | Coupled `Vm` and `phiE` | Allocates and owns `phiE` field; registered as a full factory entry |

`singleCellSolver` is registered in the parent `electroModel` table and bypasses
the myocardium-domain factory. The canonical 3D eikonal workflow selects
`myocardiumSolver eikonalSolver`, enters `electrophysiologyModel`, and builds
`eikonalMyocardiumDomain`; the legacy `myocardiumModels/eikonalSolver` class is
not a runtime-selected top-level solver.

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
| `torsoECG` | `torsoECG` | Cell-centre sampling of the unified `phiE` at electrode positions | Requires a configured `extracellularPotentialDomain`; state-provider routing handled by `electrophysicsSystemBuilder::configureECGDomains` |

---

## `conductionSystemModels/`

Concrete implementations of `conductionSystemSolver`.

| Class | Type name | Method |
|---|---|---|
| `monodomain1DSolver` | `monodomain1DSolver` | Implicit backward-Euler cable equation + ionic ODE [default] |
| `eikonalSolver1D` | `eikonalSolver1D` | Eikonal fast-marching on graph — activation times only; single param `c0` [m/s] |

**Cable equation (per edge):**

```

Cm * dVm/dt + Iion = G * d²Vm/dx²  +  Istim
where G = conductance/length

```

`monodomain1DSolver` uses the **Hines tree-elimination algorithm** for O(n) implicit solution: forward elimination leaf → root, then back-substitution root → leaf. Requires tree topology (enforced at graph load).

`eikonalSolver1D` uses a single BFS pass: `activationTime[child] = activationTime[parent] + edgeLength / c0`.

---

## `electroCouplers/`

Transfers state between domains at each timestep. Runs between domain advances in the `electrophysicsSystem` step sequence.

| Class | Role |
|---|---|
| `electroDomainCouplingEndpoints.H` | Mix-in interfaces: `tissueCouplingEndpoint` (implemented by `myocardiumDomain`) and `networkCouplingEndpoint` (implemented by conduction system domains). Provide typed access to injection targets and activation sources. |
| `electroDomainCoupler.H/C` | Base class for all couplers. Named pair of domain references with `prepareSecondaryCoupling()`, `preparePrimaryCoupling()`, `preparePostPrimaryCoupling()` hooks. |
| `pvjCoupler/pvjMapper.H/C` | Purkinje–Ventricular Junction topology mapper. Builds the spatial map between conduction-system terminal nodes and the nearest myocardium cells. |
| `pvjCoupler/pvjCoupler.H/C` | PVJ coupling-family base. Owns the shared PVJ mapper, coupling-mode parsing, and network endpoint binding. |
| `pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.H/C` | PVJ coupling with 1D-to-3D resistance model. Reads terminal `Vm`, converts it to volumetric current, and injects it into `myocardiumDomain::sourceField_`. |
| `pvjCoupler/eikonal/eikonalPvjCoupler.H/C` | PVJ coupling for activation-time models. Transfers Purkinje terminal activation times into the myocardium eikonal domain. |

Bath/extracellular-potential coupling is **not** a coupler class in this
codebase. `extracellularPotentialDomain` (an `electroStateDomain` under
`electroDomains/extracellularPotentialDomain/`) owns the global `phiE` solve
and binds a restricted view of it back into the bidomain myocardium solver
directly, replacing the older `bathDomain` / `bathECGSolver` /
`bidomainBathECGSolver` / `heartBathInterfaceCoupler` triad.

**PVJ coupling equation** (`reactionDiffusionPvjCoupler`):

```

I_pvj = (Vm_1D − Vm_3D) / R_pvj          [A/m²]
tissue source += I_pvj  (volumetric, scattered over cells within pvjRadius)
network source -= I_pvj  (bidirectional mode only)

```

`pvjMapper` spatial algorithm: for each PVJ location, find all 3D cells within `pvjRadius`, gather tissue `Vm` as a volume-weighted average, then scatter coupling current back to those cells.

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
            purkinjeGraphModelCoeffs
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
│  3. myocardiumDomain.advance(t0, dt)                          │
│       ├─ ionicModel.computeIonicCurrent()  [per cell, ODE]    │
│       ├─ myocardiumSolver.solveDiffusion*() [Vm PDE, FVM]     │
│       └─ phiE solve (bidomain only)                           │
│                                                               │
│  4. ecgDomain.advance(t0, dt)                                 │
│       └─ ecgSolver.solve()  [reads Vm from myocardium]        │
│                                                               │
│  5. write() / end() on all domains and couplers               │
└───────────────────────────────────────────────────────────────┘

```

---

## Domain field ownership

| Domain | Fields owned |
|---|---|
| `myocardiumDomain` | `Vm_`, `Iion_`, `activationTime_`, `sourceField_`, `phiE_` (bidomain only) |
| `conductionSystemDomain` | `Vm1D_`, `Iion1D_`, `activationTime_` (per node), `terminalCurrent_`, `terminalSource_` |
| `eikonalMyocardiumDomain` | `activationTime_` only — no ionic state |
| `ecgDomain` | electrode config, ECG output — reads Vm from myocardium via `electroStateProvider` |

---

## Inheritance summary

```

electroDomainInterface
    ├── myocardiumDomainInterface + electroStateProvider + tissueCouplingEndpoint
    │     ├── myocardiumDomain
    │     └── eikonalMyocardiumDomain
    ├── ecgDomain
    └── conductionSystemDomain + networkCouplingEndpoint

myocardiumSolver
    ├── monodomainSolver
    └── bidomainSolver        (owns phiE field)

electroStateProvider
    ← implemented by: myocardiumDomain
    ← consumed by:    ecgDomain, electrophysicsSystemBuilder

ElectromechanicalSignalProvider  (from couplingModels/)
    ← implemented by: ionicModel
    ← exposed by:     electroModel::provider() → myocardium → ionicModel

```
