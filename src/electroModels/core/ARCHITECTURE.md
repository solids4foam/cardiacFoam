# electroModels/core — Architecture

The `core/` directory is the orchestration layer for the multi-domain electrophysiology system.
It does not solve myocardium, Purkinje, or ECG physics directly. It does four things:

1. Selects the top-level electro model
2. Builds the assembled domain system
3. Selects the timestep orchestration scheme
4. Runs the domains and couplers in the correct order

Internally `core/` is organised into three subdirectories (`system/`,
`advanceSchemes/`, `electrophysiologyModel/`) plus flat pure-abstract interface
headers that carry no peer dependencies.

---

## Interface contracts

### `electroDomainInterface.H`

Abstract lifecycle contract that all electro domains implement. Stored as
`PtrList<electroDomainInterface>` in `electrophysicsSystem`.

| Method | Behaviour |
|---|---|
| `const Time& time() const` | Shared simulation clock |
| `void prepareTimeStep(t0, dt)` | Optional pre-advance setup (default no-op) |
| `void advance(t0, dt)` | Execute one timestep — pure virtual |
| `void write()` | Optional output and post-processing (default no-op) |
| `void end()` | Optional teardown (default no-op) |

Domain-agnostic: works with volumetric meshes, 1D graphs, or abstract solvers.

```

electroDomainInterface
├─ MyocardiumDomain  (myocardiumDomain/)
├─ ConductionSystemDomain  (conductionSystemDomain/)
└─ ECGDomain  (ecgDomain/)

```

---

### `electroVolumeFieldDomain.H`

Specialised contract for 3D FVM-based domains. A 1D graph solver implements
`electroDomainInterface` but **not** this interface.

| Method | Purpose |
|---|---|
| `const fvMesh& mesh()` | Active finite-volume mesh |
| `volScalarField& VmRef()` | Mutable transmembrane potential |
| `const volScalarField& Vm()` | Immutable Vm |
| `const volScalarField& Iion()` | Ionic current (read-only) |
| `volScalarField& sourceField()` | Coupling source term (write) |
| `const dimensionedScalar& chi()` | Surface-to-volume ratio [1/m] |
| `const dimensionedScalar& Cm()` | Membrane capacitance [F/m²] |

```

electroVolumeFieldDomain
├─ myocardiumSolver (abstract)
│  ├─ monodomainSolver
│  └─ bidomainSolver
└─ myocardiumDomain (concrete)

```

---

### `electroStateProvider.H`

Read-only, pointer-based state interface for one-way downstream coupling
(e.g., ECG reads Vm without calling back into the myocardium).

| Method | Returns |
|---|---|
| `VmPtr()` | `const volScalarField*` or `nullptr` |
| `phiEPtr()` | `const volScalarField*` or `nullptr` (monodomain → nullptr) |
| `conductivityPtr()` | `const volTensorField*` or `nullptr` |
| `intracellularConductivityPtr()` | Falls back to `conductivityPtr()` |
| `extracellularConductivityPtr()` | Falls back to `conductivityPtr()` |
| `mesh()` / `baseMesh()` | Active mesh / base mesh |
| `subsetFaceMapPtr()` / `subsetCellMapPtr()` | Non-null only for submesh strategy |

Consumer pattern:

```cpp
const auto* VmPtr = stateProvider.VmPtr();
if (VmPtr) { compute_ecg_field(*VmPtr); }

```

Breaks cyclic dependencies: ECG reads state without a callback to myocardium.

---

## Top-level selection

`electroModel::New(...)` reads `myocardiumSolver <type>` from
`constant/electroProperties` and dispatches to a registered `electroModel` subtype.

`electrophysiologyModel` is registered under `monodomainSolver`, `bidomainSolver`,
and `eikonalSolver` — the same top-level wrapper is used for all
myocardium-centred spatial workflows.

---

## Top-level model files

### `electroModel.H/C`

OpenFOAM `physicsModel` subclass. The public entry point to the entire
electrophysiology system.

- Inherits: `physicsModel`, `IOdictionary`, `electroStateProvider`

- Reads `constant/electroProperties`

- Owns the assembled `electrophysicsSystem`

- Exposes `evolve(timeValue, deltaT)` called by the main solver each timestep

- Collects per-phase performance timings

---

### `system/electrophysicsSystem.H/.C`

Domain container and advance coordinator.

- `autoPtr<myocardiumDomainInterface> myocardium_` — primary domain

- `autoPtr<electrophysicsAdvanceScheme> advanceScheme_`

- `PtrList<electroDomainInterface> conductionDomains_`

- `PtrList<ElectroDomainCoupler> conductionCouplingModels_`

- `PtrList<electroDomainInterface> ecgDomains_`

- `PtrList<ElectroDomainCoupler> ecgCouplingModels_`

The actual timestep sequence is delegated to `advanceScheme_`; this container
just holds the assembled pieces.

---

### `system/electrophysicsSystemBuilder.H/C` (namespace)

Dictionary-driven factory. Replaces hard-coded solver branching with
runtime configuration:

```cpp
// Instead of: if (type == "monodomain") { new monodomainSolver... }
auto myocardium = myocardiumDomainInterface::New(mesh, electroProperties);
system.setMyocardium(myocardium);  // works with any registered type

```

Builder functions:

| Function | Responsibility |
|---|---|
| `configureMyocardiumDomain(...)` | Instantiate myocardium domain via factory |
| `configureAdvanceScheme(...)` | Select staggered vs pimpleStaggered |
| `configureConductionDomains(...)` | Load Purkinje graph(s) and instantiate domains |
| `configureConductionCouplings(...)` | Instantiate PVJ couplers |
| `configureECGDomains(...)` | Instantiate ECG solver(s) |
| `configureECGCouplings(...)` | Instantiate myocardium-ECG couplers |

Current rules:

- `conductionNetworkDomains` and `domainCouplings` are optional

- Every conduction coupling must explicitly declare `conductionNetworkDomain <name>`

- The old single-domain fallback and `primaryDomain` field were removed

---

## Domain assembly

### Conduction (two-pass)

1. Build every entry in `conductionNetworkDomains`
2. Build every coupling in `domainCouplings`, resolving `conductionNetworkDomain <name>`
   against the already-built domain map

The two-pass exists because couplings and conduction domains are stored in
separate dictionary containers.

### ECG (single-pass)

Build every entry in `ecgDomains`. ECG domains consume myocardium state through
`electroStateProvider`; no active ECG coupling family is assembled in `core`.

---

## Advance schemes

`electrophysicsAdvanceScheme` is a runtime-selected orchestration strategy.

**Staged order (both schemes):**

1. Prepare myocardium timestep
2. Prepare conduction couplings (`prepareSecondaryCoupling`)
3. Advance conduction domains
4. Prepare myocardium couplings (`preparePrimaryCoupling`)
5. Advance myocardium
6. Prepare ECG couplings (`preparePostPrimaryCoupling`)
7. Advance ECG domains

### `advanceSchemes/staggeredElectrophysicsAdvanceScheme`

Single-pass weak coupling. Each domain sees the state from the previous
timestep. Suitable for unidirectional Purkinje → myocardium workflows.

```cpp
system.prepareConductionCouplings(t0, dt);
system.advanceConductionDomains(t0, dt);
system.prepareMyocardiumCouplings(t0, dt);
myocardium.advance(t0, dt, pimplePtr);
system.prepareECGCouplings(t0, dt);
system.advanceECGDomains(t0, dt);

```

### `advanceSchemes/pimpleStaggeredElectrophysicsAdvanceScheme`

Iterative strong coupling via OpenFOAM's PIMPLE corrector loop. The
conduction/myocardium block repeats until convergence before ECG advances.
Requires `solutionAlgorithm implicit` in electroProperties.
Suitable for bidirectional Purkinje ↔ myocardium exchange.

```cpp
while (pimplePtr->loop()) {
    system.prepareConductionCouplings(t0, dt);
    system.advanceConductionDomains(t0, dt);
    system.prepareMyocardiumCouplings(t0, dt);
    myocardium.advance(t0, dt, pimplePtr);
    system.prepareECGCouplings(t0, dt);
    system.advanceECGDomains(t0, dt);
}

```

---

## Helper files

### `overrideTypeName.H`

Drop-in replacement for OpenFOAM's `TypeName()` macro that adds `override` to
silence `-Winconsistent-missing-override` warnings in derived classes.

```cpp
// OpenFOAM standard (produces warning):
TypeName("myDerivedClass");

// This project:
OverrideTypeName("myDerivedClass");
// expands to: virtual const word& type() const override { return typeName; }

```

Used by all polymorphic solver classes in `electroModels`.

---

### `dimVoltage.H`

Shared `dimensionSet` for voltage fields so all Vm declarations are
dimensionally consistent.

```cpp
extern const dimensionSet dimVoltage;
// = [1 2 -3 0 0 -1 0]  (Volts, SI)
// = dimMass * dimArea / (pow3(dimTime) * dimCurrent)

```

---

## Coupling configuration examples

### Reaction-diffusion PVJ

```cpp
domainCouplings
{
    coupling1
    {
        electroDomainCoupler    reactionDiffusionPvjCoupler;
        conductionNetworkDomain purkinjeNetwork;
        couplingMode            unidirectional;
        pvjRadius               6e-4;
        rPvj                    500.0;
    }
}

```

### Eikonal PVJ

```cpp
domainCouplings
{
    coupling1
    {
        electroDomainCoupler    eikonalPvjCoupler;
        conductionNetworkDomain purkinjeNetwork;
        couplingMode            unidirectional;
        pvjRadius               6e-4;
    }
}

```

`bidirectional` is parsed but explicitly rejected as in development for
`eikonalPvjCoupler`.

---

## Extension guidance

The intended rule when extending `core`:

- `core` owns orchestration only

- domain layers own domain-family selection and state

- solver folders own numerical kernels

- coupler folders own exchange laws

New solver-family branching should go into a domain-layer factory, not into the
builder.

---

## Architecture diagram

```

┌─────────────────────────────────────────────────────────────┐
│ electroModel (top-level facade)                             │
│  • Reads constant/electroProperties                         │
│  • Owns electrophysicsSystem                                │
│  • Implements electroStateProvider (exposes Vm, phiE, σ)   │
│  • Entry point: evolve(t, dt)                               │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ electrophysicsSystem (domain container)                     │
│  • myocardium (implements electroVolumeFieldDomain)         │
│  • advanceScheme (staggered or pimpleStaggered)             │
│  • conductionDomains (Purkinje graphs)                      │
│  • ecgDomains (ECG solvers)                                 │
│  • couplers (Purkinje↔myocardium, myocardium→ECG)          │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌──────────────────────────────┬──────────────────────────────┐
│ electrophysicsAdvanceScheme  │ electrophysicsSystemBuilder  │
│ (abstract strategy)          │ (dictionary-driven factory)  │
│ • advance()                  │ • configureMyocardium()      │
├──────────────────────────────┤ • configureConduction()      │
│ Implementations:             │ • configureECG()             │
│ • staggered (weak)           │ • configureCouplers()        │
│ • pimpleStaggered (strong)   │                              │
└──────────────────────────────┴──────────────────────────────┘

```

## File dependency overview

```

electroModel.H
  ├─ electroStateProvider.H
  ├─ system/electrophysicsSystem.H
  │  ├─ electroDomainInterface.H
  │  │  ├─ myocardiumDomain (implements)
  │  │  ├─ ConductionSystemDomain (implements)
  │  │  └─ ECGDomain (implements)
  │  ├─ advanceSchemes/electrophysicsAdvanceScheme.H
  │  │  ├─ advanceSchemes/staggered/staggeredElectrophysicsAdvanceScheme.H
  │  │  └─ advanceSchemes/pimpleStaggered/pimpleStaggeredElectrophysicsAdvanceScheme.H
  │  └─ ElectroDomainCoupler.H
  └─ system/electrophysicsSystemBuilder.H
     ├─ electroVolumeFieldDomain.H
     └─ dimVoltage.H

overrideTypeName.H  (used by all polymorphic solver classes)

```

---

## Key design principles

1. **Domain separation** — each domain implements `electroDomainInterface` and advances independently.
2. **Dependency inversion** — domains depend on abstract interfaces (`electroStateProvider`, `electroDomainInterface`), not concrete solvers.
3. **One-way coupling** — post-domains (ECG) read state via `electroStateProvider` pointers; they never call back into the primary domain.
4. **Runtime selection** — `electrophysicsSystemBuilder` enables dictionary-driven instantiation instead of hard-coded branching.
5. **Pluggable orchestration** — `electrophysicsAdvanceScheme` subclasses allow swapping weak vs strong coupling without changing any domain code.
6. **Type safety** — `electroVolumeFieldDomain` enforces that 3D solvers provide volumetric fields; graph solvers cannot accidentally fill a 3D-mesh role.
