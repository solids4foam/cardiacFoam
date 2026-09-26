# electroModels/core — Architecture

The `core/` directory is the orchestration layer for the multi-domain electrophysiology system.
It does not solve myocardium, Purkinje, or ECG physics directly. It does four things:

1. Selects the top-level electro model
2. Builds the assembled domain system
3. Selects the timestep orchestration scheme
4. Runs the domains and couplers in the correct order

Internally `core/` is organised into four subdirectories (`system/`,
`advanceSchemes/`, `electrophysiologyModel/`, `verificationModels/`) plus flat
pure-abstract interface headers that carry no peer dependencies.
`verificationModels/` holds the abstract verifier base classes
(`electroVerificationModel`, `couplingVerificationModel`, `ecgVerificationModel`,
`eikonalVerificationModel`, `graphVerificationModel`) that the concrete
verifiers in `src/verificationModels` inherit from.

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
├─ myocardiumDomain  (myocardiumDomain/)
├─ conductionSystemDomain  (conductionSystemDomain/)
└─ ecgDomain  (ecgDomain/)

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
| `const volScalarField* IionOldPtr()` / `IionOldOldPtr()` | Previous timestep(s)' ionic current, for solvers that extrapolate it (default: `nullptr`) |
| `volScalarField& sourceField()` | Coupling source term (write) |
| `volScalarField* implicitSourceCoeffPtr()` (const and non-const) | Optional implicit source coefficient (default: `nullptr`) |
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
| `activationTimePtr()` | `const volScalarField*` or `nullptr` |
| `phiEPtr()` | `const volScalarField*` or `nullptr` (monodomain → nullptr) |
| `conductivityPtr()` | `const volTensorField*` or `nullptr` |
| `chiPtr()` / `CmPtr()` / `c0Ptr()` | `const dimensionedScalar*` or `nullptr` |
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

### `electroStateDomain.H`

Combined contract for domains that both advance in time and expose read-only
state: `public electroDomainInterface, public electroStateProvider`. Used for
the bath/extracellular potential domain, held by `electrophysicsSystem` as
`autoPtr<electroStateDomain> potentialDomain_` and built by
`configureBathPotentialDomain(...)`.

```

electroStateDomain
├─ electroDomainInterface (implements)
├─ electroStateProvider (implements)
└─ extracellularPotentialDomain (concrete)

```

---

## Top-level selection

`electroModel::New(...)` reads `myocardiumSolver <type>` from
`constant/electroProperties` and dispatches to a registered `electroModel` subtype.

`electrophysiologyModel` is registered under `monodomainSolver`, `bidomainSolver`,
and `eikonalSolver` — the same top-level wrapper is used for all
myocardium-centred spatial workflows.

Inside `electrophysiologyModel`, `myocardiumDomainInterface::New(...)` performs
the secondary myocardium-domain dispatch. `monodomainSolver` and
`bidomainSolver` resolve through the `myocardiumSolver` table; `eikonalSolver`
builds `eikonalMyocardiumDomain`. `singleCellSolver` is registered directly in
the parent `electroModel` table and does not enter this secondary dispatch.

---

## Top-level model files

### `electroModel.H/C`

OpenFOAM `physicsModel` subclass. The public entry point to the entire
electrophysiology system.

- Inherits: `physicsModel`, `IOdictionary`, `electroStateProvider`

- Reads `constant/electroProperties`

- Owns the assembled `electrophysicsSystem`

- Exposes `evolve()` called by the main solver each timestep (the `physicsModel`
  interface takes no arguments; timing comes from the shared `Time` object)

- Collects per-phase performance timings

---

### `electrophysiologyModel/printElectrophysiologySummary.H`

Startup banner and configuration summary parser.
Included in the `electrophysiologyModel` constructor to echo the selected solver hierarchy, advance scheme, active ionic models, tissue heterogeneities, and spatial dimensions to the terminal prior to execution.

---

### `system/electrophysicsSystem.H/.C`

Domain container and advance coordinator.

- `autoPtr<myocardiumDomainInterface> myocardium_` — primary domain

- `autoPtr<electroStateDomain> potentialDomain_` — optional bath/extracellular
  potential domain, only set for unified-phiE (bath-coupled bidomain) cases

- `autoPtr<electrophysicsAdvanceScheme> advanceScheme_`

- `PtrList<electroDomainInterface> conductionDomains_`

- `PtrList<electroDomainCoupler> conductionCouplingModels_`

- `PtrList<electroDomainInterface> ecgDomains_`

- `PtrList<electroDomainCoupler> ecgCouplingModels_`

The actual timestep sequence is delegated to `advanceScheme_`; this container
just holds the assembled pieces.

---

### `system/electrophysicsSystemBuilder.H/C` (namespace)

Dictionary-driven factory.

```cpp
auto myocardium = myocardiumDomainInterface::New(mesh, electroProperties);
system.setMyocardium(myocardium);
```

Builder functions, called in this order by the caller (`electrophysiologyModel`'s
constructor) — the order matters, since later calls depend on state set up by
earlier ones:

| Function | Responsibility |
|---|---|
| `configureMyocardiumDomain(...)` | Instantiate myocardium domain via factory |
| `configureAdvanceScheme(...)` | Select the runtime advance scheme |
| `configureBathPotentialDomain(...)` | Optionally instantiate a standalone bath/extracellular potential domain (an `electroStateDomain`, concretely `extracellularPotentialDomain`); requires the myocardium domain to already be configured |
| `configureConductionDomains(...)` | Load Purkinje graph(s), instantiate domains, and build their PVJ couplings in the same pass |
| `configureECGDomains(...)` | Instantiate ECG domain(s), reading myocardium (and optionally bath) state via `electroStateProvider`, and build any per-domain ECG couplings in the same pass |

No `configure*Couplings(...)` functions exist.

`configureBathPotentialDomain(...)` instantiates an optional standalone
domain, gated on whether `bathPotentialDomain` is present in
`electroProperties`. It builds no couplings. It requires the myocardium
domain to already be configured (`system.hasMyocardium()`). It runs before
`configureECGDomains(...)`: the potential domain it builds is one of the
state providers `configureECGDomains(...)` routes a `torsoECG` domain to.

Coupling instantiation for conduction and for ECG is built inline, in two
different structures:

- conduction (PVJ) couplings: built inside `configureConductionDomains(...)`
  from a top-level `domainCouplings` block; each entry resolves against the
  conduction-domain map via `conductionNetworkDomain <name>`
- ECG couplings: built inside `configureECGDomains(...)` from an optional
  `coupling` subdict nested under each `ecgDomains.<name>` entry; there is
  no top-level `domainCouplings`-equivalent container for ECG

Current rules:

- `conductionNetworkDomains` and `domainCouplings` are optional
- Every conduction coupling must explicitly declare `conductionNetworkDomain <name>`
- `bathPotentialDomain` is optional; when present it is read from
  `electroProperties.subDict("bathPotentialDomain")`
- Each `ecgDomains.<name>` entry's `coupling` subdict is optional; when
  present, `configureECGDomains(...)` requires a myocardium domain to exist

---

## Domain assembly

### Bath potential domain (single, optional)

Build the one `bathPotentialDomain` entry if present. Independent of the
conduction and ECG passes; only depends on the myocardium domain already
being built. Built before conduction, since a `torsoECG`-class ECG domain
built afterwards may need it.

### Conduction (two-pass)

1. Build every entry in `conductionNetworkDomains`
2. Build every coupling in `domainCouplings`, resolving `conductionNetworkDomain <name>`
   against the already-built domain map

The two-pass exists because couplings and conduction domains are stored in
separate dictionary containers.

### ECG (two-pass)

1. Build every entry in `ecgDomains`. Each domain's `ecgSolver` (`pseudoECG`,
   `eikonalECG`, or `torsoECG`) determines which `electroStateProvider` it
   reads from — `torsoECG` requires the bath potential domain from the
   previous step.
2. For each ECG domain whose dict declares a `coupling` subdict, build a
   coupler (`electroDomainCoupler::New(myocardium, ecgDomain, couplingDict)`)
   and append it to `ecgCouplingModels_`. Domains without a `coupling`
   subdict get none.

ECG couplings are declared per-domain, under each `ecgDomains.<name>`
entry's `coupling` subdict, rather than in a separate top-level container
like conduction's `domainCouplings`.

---

## Advance schemes

`electrophysicsAdvanceScheme` is a runtime-selected orchestration strategy.

**Staged order** (method names as declared on `electrophysicsSystem`):

1. Prepare myocardium timestep (`myocardium.prepareTimeStep(...)`)
2. Prepare conduction couplings (`prepareConductionCouplings(...)`)
3. Advance conduction domains (`advanceConductionDomains(...)`)
4. Prepare myocardium couplings (`prepareMyocardiumCouplings(...)`)
5. Advance the primary domain — either:
   - no potential domain: `myocardium.advance(t0, dt, pimplePtr)`, or
   - a potential domain is set (`hasPotentialDomain()`): a split
     reaction/diffusion sequence coupling `myocardium` and
     `potentialDomain_` (see below)
6. Prepare ECG couplings (`prepareECGCouplings(...)`)
7. Advance ECG domains (`advanceECGDomains(...)`)

### `advanceSchemes/staggeredElectrophysicsAdvanceScheme`

Single-pass weak coupling. Each domain sees the state from the previous
timestep. Suitable for unidirectional Purkinje → myocardium workflows.

```cpp
myocardium.prepareTimeStep(t0, dt);
system.prepareConductionCouplings(t0, dt);
system.advanceConductionDomains(t0, dt);
system.prepareMyocardiumCouplings(t0, dt);
myocardium.advance(t0, dt, pimplePtr);
system.prepareECGCouplings(t0, dt);
system.advanceECGDomains(t0, dt);

```

**Bath/unified-phiE branch:** when `system.hasPotentialDomain()` (a bath
potential domain was configured via `configureBathPotentialDomain(...)`),
step 5 above is replaced by a split reaction/diffusion sequence instead of
the plain `myocardium.advance(...)` call — this requires a myocardium
domain that `supportsSplitReactionDiffusion()`:

```cpp
myocardium.solveReactionStep(t0, dt);
system.preparePotentialDomain(t0, dt);

if (bathPredictorCorrector_)   // default true
{
    // Predict Vm, update phiE, then correct Vm.
    myocardium.solveDiffusionStepOnce(t0, dt, pimplePtr);
    system.advancePotentialDomain(t0, dt);
    myocardium.solveDiffusionStepOnce(t0, dt, pimplePtr);
}
else
{
    // Update phiE, then solve Vm with the updated phiE held fixed.
    system.advancePotentialDomain(t0, dt);
    myocardium.solveDiffusionStep(t0, dt, pimplePtr);
}

myocardium.finalizeDiffusionStep();

```

This is how the bath/extracellular-potential coupling (unified phiE) is
implemented: `staggeredElectrophysicsAdvanceScheme` owns a
`bathPredictorCorrector_` switch (dictionary key
`bathPredictorCorrector`, default `true`) that picks between the two
sub-variants above.

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
        pvjCouplingScheme       implicit; // optional: explicit | implicit
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

`bidirectional` is implemented for `eikonalPvjCoupler`: it gathers myocardial
activation times onto the network terminals via `mapper_.gatherActivationTimes`.

---

## Extension rules

When extending `core`:

- `core` owns orchestration only
- Domain layers own domain-family selection and state
- Solver folders own numerical kernels
- Coupler folders own exchange laws

New solver-family branching belongs in a domain-layer factory, not in the builder.

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
│  • potentialDomain (optional, implements electroStateDomain)│
│  • advanceScheme                                             │
│  • conductionDomains (Purkinje graphs)                       │
│  • ecgDomains (ECG solvers)                                  │
│  • couplers (Purkinje↔myocardium, myocardium→ECG)          │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌──────────────────────────────┬──────────────────────────────┐
│ electrophysicsAdvanceScheme  │ electrophysicsSystemBuilder  │
│ (abstract strategy)          │ (dictionary-driven factory)  │
│ • advance()                  │ • configureMyocardiumDomain()│
├──────────────────────────────┤ • configureAdvanceScheme()   │
│ Implementations:             │ • configureConductionDomains()│
│ • staggered (weak, with      │ • configureBathPotentialDomain()│
│   optional bath coupling)    │ • configureECGDomains()      │
└──────────────────────────────┴──────────────────────────────┘

```

## File dependency overview

```

electroModel.H
  ├─ electromechanicalSignalProvider.H
  ├─ electroStateProvider.H
  └─ system/electrophysicsSystem.H
     ├─ electroDomainInterface.H
     │  ├─ myocardiumDomain (implements)
     │  ├─ conductionSystemDomain (implements)
     │  └─ ecgDomain (implements)
     ├─ electroStateDomain.H  (electroDomainInterface + electroStateProvider)
     │  └─ extracellularPotentialDomain (implements; the bath/potential domain)
     ├─ electroDomainCoupler.H
     ├─ myocardiumDomainInterface.H
     └─ advanceSchemes/electrophysicsAdvanceScheme.H
        └─ advanceSchemes/staggered/staggeredElectrophysicsAdvanceScheme.H

system/electrophysicsSystemBuilder.H
  ├─ system/electrophysicsSystem.H  (as above)
  └─ electroStateProvider.H

overrideTypeName.H  (used by all polymorphic solver classes)
dimVoltage.H / electroVolumeFieldDomain.H  (used by the myocardium domain
  layer, e.g. myocardiumDomain.H / myocardiumSolver.H — not included from
  core's builder or system files directly)

```

---

## Key design principles

1. **Domain separation** — each domain implements `electroDomainInterface` and advances independently.
2. **Dependency inversion** — domains depend on abstract interfaces (`electroStateProvider`, `electroDomainInterface`), not concrete solvers.
3. **One-way coupling** — post-domains (ECG) read state via `electroStateProvider` pointers; they never call back into the primary domain.
4. **Runtime selection** — `electrophysicsSystemBuilder` enables dictionary-driven instantiation instead of hard-coded branching.
5. **Pluggable orchestration** — `electrophysicsAdvanceScheme` subclasses allow swapping weak vs strong coupling without changing any domain code.
6. **Type safety** — `electroVolumeFieldDomain` enforces that 3D solvers provide volumetric fields; graph solvers cannot accidentally fill a 3D-mesh role.
