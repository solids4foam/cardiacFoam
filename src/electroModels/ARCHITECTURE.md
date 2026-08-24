# electroModels

Multi-domain spatial electrophysiology framework for OpenFOAM-based cardiac simulation.

Compiled as `libelectroModels`. Provides co-simulation infrastructure for three cardiac subsystems — myocardium, body-surface ECG, and the conduction system — each on a separate mesh, coupled through typed interfaces at every timestep.

## Runtime flow and public selectors

```text
physicsModel::New()             physicsProperties: type electroModel
  -> electroModel::New()        electroProperties: myocardiumSolver <name>
     -> electrophysiologyModel  monodomainSolver | bidomainSolver | eikonalSolver
        -> myocardium domain
           -> myocardiumSolver  monodomainSolver | bidomainSolver
```

The eikonal spatial entry uses `eikonalMyocardiumDomain` rather than the
`myocardiumSolver` kernel table. `singleCellSolver` registers directly in the
parent `electroModel` table and bypasses this spatial assembly. For spatial
entries, `<name>Coeffs` is the builder's configuration root; optional public
blocks include `conductionNetworkDomains`, `ecgDomains`,
`bathPotentialDomain`, and `domainCouplings`. Domain-specific selectors include
`conductionSystemSolver` and `ecgSolver`.

This complete electrophysiology path is built in both full solids4foam and
lightweight modes. The build-mode boundary is outside this library:
`electroMechanicalModels` is built only in full mode.

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

The framework backbone. Owns orchestration only; domain-family selection is delegated to domain-layer factories.

| File | Description |
|---|---|
| `electroModel.H/C` | OpenFOAM `physicsModel` subclass. Solver-level entry point. Owns the `electrophysicsSystem`, drives the time loop via `evolve()`, and delegates all stepping to the system. |
| `system/electrophysicsSystem.H/C` | Holds the assembled domain hierarchy: typed `myocardium_`, `conductionDomains_`, `ecgDomains_`, and coupling model lists. Implements the ordered advance → couple → write sequence. |
| `system/electrophysicsSystemBuilder.H/C` | Builds the full system from `electroProperties`. Instantiates the myocardium domain through `myocardiumDomainInterface::New(...)`, then assembles conduction and ECG-stage domains. |
| `electrophysiologyModel/` | Top-level myocardium-centered orchestration wrapper registered under `monodomainSolver`, `bidomainSolver`, and `eikonalSolver`. Owns ionic model, verification model, and field export lists. |
| `advanceSchemes/electrophysicsAdvanceScheme.H/C` | Abstract time-advance strategy. Defines when the ionic solve and the diffusion solve are called relative to each other. |
| `advanceSchemes/staggered/` | Staggered (operator-split) time advance: single-pass weak coupling. |
| `electroDomainInterface.H` | Minimal lifecycle contract that all electro domains implement. Pure virtual: `time()`, `advance(t0, dt)`. Optional no-ops: `prepareTimeStep()`, `write()`, `end()`. |
| `electroStateProvider.H` | Read-only field interface implemented by domains that expose fields upstream: `VmPtr()`, `phiEPtr()`, `conductivityPtr()`. Consumed by ECG solver and the builder. |
| `dimVoltage.H` | Shared dimension set for voltage fields (`[1 2 -3 0 0 -1 0]`, i.e. V). |
| `overrideTypeName.H` | Macro to assign a lowercase runtime type name, used with `addToRunTimeSelectionTable`. |

---

## `electroDomains/`

Domains own specific meshes and fields, communicating via `electroStateProvider` and `electroCouplers`.

### `myocardiumDomain/`

Primary 3D cardiac region. Selects concrete solvers via `myocardiumDomainInterface`.

**`myocardiumDomain`**

- Inherits: `electroDomainInterface`, `tissueCouplingEndpoint`, `electroStateProvider`
- Owns: `autoPtr<myocardiumSolver>`, `ionicModel&`

- Fields: `Vm_` (transmembrane voltage), `sourceField_` (external current injection)

- On construction: detects whether the solver provides `phiE` (bidomain case) and binds it

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

Concrete reaction-diffusion implementations of `myocardiumSolver`.

| Class | Type name | PDE / method | Notes |
|---|---|---|---|
| `monodomainSolver` | `monodomainSolver` | `∂Vm/∂t − ∇·(σᵢ∇Vm) = Iion` | Standard single-domain FVM |
| `bidomainSolver` | `bidomainSolver` | Coupled `Vm` and `phiE` | Owns a local `phiE`; with a bath domain, synchronizes it from the global `phiE` owner |

`singleCellSolver` is registered in the parent `electroModel` table and bypasses the myocardium-domain factory. The canonical eikonal workflow selects `myocardiumSolver eikonalSolver` and builds `eikonalMyocardiumDomain`.

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
| `eikonalECG` | `eikonalECG` | Template-voltage surrogate: reconstructs `Vm(x,t) = U(t - psi(x))` from endo/mid/epi single-cell templates, then pseudo-ECG integrates | For eikonal activation-time workflows; no reaction-diffusion solve per step. See [ecgModels/eikonalECG/README.md](ecgModels/eikonalECG/README.md) |

---

## `conductionSystemModels/`

Concrete implementations of `conductionSystemSolver`.

| Class | Type name | Method |
|---|---|---|
| `monodomain1DSolver` | `monodomain1DSolver` | Implicit backward-Euler cable equation + ionic ODE [default] |
| `eikonalSolver1D` | `eikonalSolver1D` | Eikonal fast-marching on graph — activation times only; single param `c0` [m/s] |
| `restitutionEikonalSolver1D` | `restitutionEikonalSolver1D` | Re-excitable activation solver. Beat-to-beat interval logic per node with CV(DI) restitution and constant `apdNominal`; reports block, wavebreak, short-DI and minimum-DI diagnostics. See [conductionSystemModels/README.md](conductionSystemModels/README.md) |

**Cable equation (per edge):**

```

Cm * dVm/dt + Iion = G * d²Vm/dx²  +  Istim
where G = conductance/length

```

`monodomain1DSolver` uses the Hines tree-elimination algorithm — O(n). Requires tree topology.

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
| `pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.H/C` | PVJ coupling with 1D-to-3D resistance model. Reads terminal `Vm`, converts it to volumetric current, and injects it explicitly into `myocardiumDomain::sourceField_` or, with `pvjCouplingScheme implicit`, splits the tissue-voltage term into the myocardium Vm matrix diagonal. |
| `pvjCoupler/eikonal/eikonalPvjCoupler.H/C` | PVJ coupling for activation-time models. Transfers Purkinje terminal activation times into the myocardium eikonal domain. |
| `pvjCoupler/eikonalMonodomain/eikonalMonodomainPvjCoupler.H/C` | Eikonal-network to monodomain-tissue PVJ coupling. Anterograde transfer uses a voltage template; under `couplingMode bidirectional` it additionally gathers myocardial activation times back onto the network terminals (retrograde 3D-to-1D). This is the only coupler with a working bidirectional path. |

Bath/extracellular-potential coupling is **not** a coupler class. `extracellularPotentialDomain` (an `electroStateDomain` under `electroDomains/extracellularPotentialDomain/`) owns the global `phiE` solve and binds a restricted view of it directly into the bidomain myocardium solver.

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
                conductionSystemSolver  monodomain1DSolver;  // default; or: eikonalSolver1D | restitutionEikonalSolver1D
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

## Domain field ownership

| Domain | Fields owned |
|---|---|
| `myocardiumDomain` | `Vm_`, `Iion_`, `activationTime_`, `sourceField_`; its bidomain kernel owns local `phiE` and can synchronize it from the bath domain |
| `conductionSystemDomain` | `Vm1D_`, `Iion1D_`, `activationTime_` (per node), `terminalCurrent_`, `terminalSource_` |
| `eikonalMyocardiumDomain` | `activationTime_` only — no ionic state |
| `ecgDomain` | electrode config, ECG output — reads Vm from myocardium via `electroStateProvider` |
| `extracellularPotentialDomain` | global `phiE`, `sigmaTotal`, and scattered `VmGlobal` |

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
    └── bidomainSolver        (owns local phiE; may bind global phiE)

electroStateProvider
    ← implemented by: myocardiumDomain
    ← consumed by:    ecgDomain, electrophysicsSystemBuilder

ElectromechanicalSignalProvider  (from couplingModels/)
    ← implemented by: ionicModel
    ← exposed by:     electroModel::provider() → myocardium → ionicModel

```

## Source boundaries

This directory is hand-maintained project code and is listed explicitly in
`Make/files`. It consumes `ionicModels`, `genericWriter`, and the selected
`physicsModel` interface. Generated ionic equation headers live under
`src/ionicModels`, not here. `modules/solids4foam` is external submodule
content; compatibility adaptations belong in cardiacFoam or its owned
lightweight `modules/physicsModel` layer.
