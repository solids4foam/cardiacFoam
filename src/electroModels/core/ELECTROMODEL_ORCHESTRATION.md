# electroModel — Core Orchestration and Multi-Domain System Architecture

## 1. Overview

The `electroModel` class is the **top-level physics model entry point** for the cardiac electrophysiology solver. It inherits from OpenFOAM's `physicsModel` interface and acts as the facade that:

- Reads `constant/electroProperties` configuration dictionary
- Owns the `electrophysicsSystem` — the multi-domain container
- Drives the main time loop via the `evolve()` interface
- Delegates all physics stepping to the system's advance scheme

This document explains how the orchestration layer binds together three distinct physical domains (myocardium, conduction system, ECG) into a coherent simulation.

---

## 2. The electroModel Interface

### Constructor and Initialization

```cpp
electroModel(const fvMesh& mesh, const IOdictionary& dict);
```

On construction:
1. **Reads electroProperties** from `constant/` directory
2. **Instantiates electrophysicsSystem** via `electrophysicsSystemBuilder`
3. **Selects advance scheme** from runtime table (default: `staggeredElectrophysicsAdvanceScheme`)
4. **Configures domains** based on presence/absence of optional coupling blocks

### Main Evolution Loop

```cpp
virtual bool evolve();
```

Called by cardiacFoam's time-integration loop. Delegates to the advance scheme:

```cpp
bool success = advanceScheme_->advance(runTime_.value(), deltaT, system_, pimplePtr_, timings_);
```

---

## 3. The electrophysicsSystem Container

The **orchestration workhorse**. Holds:

### Primary Domain
```cpp
autoPtr<MyocardiumDomain> myocardium_;  // Always present, 3D tissue
```

### Upstream Domains (Pre-Primary)
```cpp
PtrList<electrophysicsSystem::upstreamDomain> upstreamDomains_;      // 0 or 1: Purkinje
PtrList<ElectroDomainCoupler> upstreamCouplingModels_;               // 0 or 1: pvjResistanceCoupler
```

These advance **before** myocardium and couple into it (e.g., Purkinje activation current → tissue stimulus).

### Downstream Domains (Post-Primary)
```cpp
PtrList<electroDomainInterface> downstreamDomains_;                 // 0 or 1: ECG
PtrList<ElectroDomainCoupler> downstreamCouplingModels_;            // Usually empty
```

These advance **after** myocardium and read its voltage field (one-way data flow).

### Time Advance Scheme
```cpp
autoPtr<electrophysicsAdvanceScheme> advanceScheme_;
```

Encapsulates the stepping strategy (weak/strong coupling, explicit/implicit).

---

## 4. The Advance Scheme Strategy Pattern

### Base Class: electrophysicsAdvanceScheme

Defines the pure virtual interface:

```cpp
virtual bool advance
(
    scalar t0,
    scalar dt,
    electrophysicsSystem& system,
    pimpleControl* pimplePtr,     // Passed from implicit solver
    electrophysicsAdvanceTimings& timings
) = 0;
```

The scheme orchestrates:
1. When domains prepare and couple
2. In what order domains step
3. How many iterations (if implicit)

### Concrete Schemes

#### staggeredElectrophysicsAdvanceScheme (Explicit Staggered)

```
Prepare upstream coupling (read Vm(n))
   ↓
Advance upstream (Purkinje to Vm(n+1))
   ↓
Prepare primary coupling (compute coupling term with new upstream voltage)
   ↓
Advance primary (myocardium to Vm(n+1) with updated stimulus)
   ↓
Prepare + Advance downstream (ECG reads finalized tissue voltage)
```

**Characteristics:**
- ✅ Computationally efficient (single pass per time step)
- ✅ Stable for unidirectional coupling (Purkinje → tissue)
- ⚠️ Weak coupling introduces lag; stiff coupling can oscillate

#### pimpleStaggeredElectrophysicsAdvanceScheme (Implicit Staggered with PIMPLE)

```
while (pimplePtr->loop())  // Iterate within time step
{
    Prepare upstream coupling
       ↓
    Advance upstream (with pimplePtr)
       ↓
    Prepare primary coupling
       ↓
    Advance primary (with pimplePtr)
}
   ↓
Prepare + Advance downstream
```

**Characteristics:**
- ✅ Strongly coupled (multiple iterations → convergence)
- ✅ Stable for bidirectional coupling (tissue → Purkinje feedback)
- ⚠️ Higher per-timestep cost, but larger stable Δt possible

---

## 5. Configuration Flow

### Dictionary Parsing (electroProperties)

```
electroProperties
├── myocardiumSolver <type>              // Selects 3D domain solver
├── <type>Coeffs { ... }                  // Domain-specific parameters
├── electrophysicsAdvanceScheme <type>   // Time-stepping strategy
├── conductionNetworkDomains { ... }      // Optional upstream domains (Purkinje)
├── domainCouplings { ... }               // Optional couplers (pvjResistanceCoupler)
└── ecgDomains { ... }                    // Optional downstream domains (ECG)
```

### Builder Pattern: electrophysicsSystemBuilder

Free functions that instantiate the system:

```cpp
std::unique_ptr<electrophysicsSystem> buildSystem
(
    const fvMesh& mesh,
    const dictionary& electroProperties
);
```

The builder:
1. Detects which optional domains are present (via dictionary keys)
2. Calls domain constructors in sequence
3. Instantiates couplers if coupling blocks exist
4. Wires everything into the system

**Example:**

```cpp
// If conductionNetworkDomains is present:
buildUpstreamDomains(mesh, dict)
    → purkinjeNetworkModel(mesh, couplingDict)
    → GraphConductionSystemDomain owns GraphConductionSystemSolver

// If domainCouplings is present:
buildCouplers(mesh, dict)
    → pvjResistanceCoupler(mesh, couplingDict)
```

---

## 6. Data Flow in One Time Step (Example: Staggered Scheme)

Assume: monodomainSolver + purkinjeNetworkModel + pvjResistanceCoupler

### Step 1: Prepare Upstream

```
pvjResistanceCoupler::prepareSecondaryCoupling(t0, dt)
    → Reads Vm(n) from myocardium (stale)
    → Reads Vp(n) from Purkinje
    → Computes IPM = (Vp - Vm) / Rpvj
```

### Step 2: Advance Upstream

```
purkinjeNetworkModel::advance(t0, dt)
    → Reads coupling current IPM from coupler
    → Solves 1D activation on graph
    → Updates Vp(n+1)
```

### Step 3: Prepare Primary

```
pvjResistanceCoupler::preparePrimaryCoupling(t0, dt)
    → Reads updated Vp(n+1) from Purkinje
    → Reads stale Vm(n) from myocardium (not yet updated)
    → Recomputes IPM = (Vp(n+1) - Vm(n)) / Rpvj
    → Deposits into myocardium sourceField
```

### Step 4: Advance Primary

```
MyocardiumDomain::advance(t0, dt)
    → Solves 3D monodomain PDE with IPM in sourceField
    → Updates Vm(n+1)
```

### Step 5: Advance Downstream

```
ecgDomain::advance(t0, dt)
    → Reads finalized Vm(n+1) from myocardium
    → Solves passive potential field (Laplace/Poisson)
    → Computes ECG at electrode sites
```

---

## 7. Key Design Decisions

### Why Not Fully Implicit?

A fully implicit (monolithic) block-coupled system would require:
- Single massive matrix containing both tissue and Purkinje
- Complex block-matrix extensions to OpenFOAM
- Significant memory overhead

**Decision:** Segregated iterative coupling via PIMPLE is a middle ground:
- Reuses OpenFOAM's standard fvMatrix machinery
- Allows independent domain solvers
- Achieves implicit stability with modest overhead

### Why Separate Advance and Prepare?

The **two-phase coupling** (prepare ↔ advance) allows:
- Couplers to query stale fields from previous sub-step
- Explicit evaluation of coupling terms
- Flexibility to implement semi-implicit extensions (diagonal matrix modification)

### Why Optional Domains?

Users may need:
- **Myocardium only** → skip Purkinje, skip ECG
- **Myocardium + Purkinje** → unidirectional driving (fast)
- **Myocardium + Purkinje (bidirectional)** → requires pimple scheme + implicit
- **Myocardium + ECG** → one-way monitoring

Dictionary-driven construction avoids code branching and allows runtime flexibility.

---

## 8. Extension Points

### Adding a New Domain Type

1. Inherit `electrophysicsSystem::upstreamDomain` or `downstreamDomain`
2. Implement `advance()`, optional `prepareTimeStep()`, `write()`, `end()`
3. Register in runtime table via `addToRunTimeSelectionTable()`
4. Add builder logic to `electrophysicsSystemBuilder`
5. Document in `electroProperties.template`

### Adding a New Advance Scheme

1. Inherit `electrophysicsAdvanceScheme`
2. Implement `advance()` with custom stepping order/iteration
3. Register in runtime table
4. Document in template and dict_entries.py

**Example:** A future **semi-implicit scheme** could:
- Compute the resistive coefficient 1/Rpvj
- Pass to both domains' matrix systems
- Treat neighbor voltage explicitly on RHS
- Single pass per time step with better stability than explicit staggered

---

## 9. Summary

| Component | Role | Flexibility |
|-----------|------|-------------|
| `electroModel` | Entry point, I/O, time-stepping control | Fixed (interface) |
| `electrophysicsSystem` | Multi-domain container, orchestration | Via dictionary |
| `electrophysicsAdvanceScheme` | Step strategy (order, iteration count) | Runtime selectable |
| Domains (myocardium, Purkinje, ECG) | Physics solvers | Runtime selectable |
| Couplers | Inter-domain data exchange | Runtime selectable |

The architecture cleanly separates **what to solve** (domains) from **how to solve** (scheme), enabling researchers to quickly experiment with different coupling approaches without modifying core orchestration logic.

