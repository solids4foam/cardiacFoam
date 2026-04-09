# src/ Documentation Guide — Complete Architecture Overview

## Overview

The `src/` directory contains all **source code libraries** that implement cardiacFoam's physics and algorithms. This guide maps you to the detailed documentation for each major module.

---

## Module Documentation Map

### 1. Electrophysiology (electroModels/)

**Solves the electrical activation of cardiac tissue.**

| Document | Scope | Covers |
|----------|-------|--------|
| **[ARCHITECTURE.md](./electroModels/ARCHITECTURE.md)** | All domains & solvers | Directory structure, domain types, solver hierarchy |
| **[README.md](./electroModels/README.md)** | High-level overview | Components, runtime selection, build process |
| **[ELECTROMODEL_ORCHESTRATION.md](./electroModels/core/ELECTROMODEL_ORCHESTRATION.md)** | Core coordination | electroModel facade, electrophysicsSystem, advance schemes |
| **[electroCouplers/README.md](./electroModels/electroCouplers/README.md)** | Domain coupling | PVJ coupling endpoints, mapper, resistance coupler, exchange modes |
| **[electroDomains/README.md](./electroModels/electroDomains/README.md)** | Domain interfaces | MyocardiumDomain, ECGDomain, ConductionSystemDomain contracts |

**Key Files to Know:**
- `core/electroModel/electroModel.H` — Entry point
- `core/schemes/` — Time-stepping strategies (staggered, pimpleStaggered)
- `electroDomains/myocardiumDomain/` — 3D tissue solver
- `electroDomains/conductionSystemDomain/` — 1D Purkinje network
- `myocardiumModels/` — 3D PDE solvers (mono/bidomain/eikonal)
- `electroCouplers/` — Domain coupling (PVJ resistance model)

---

### 2. Ionic Models (ionicModels/)

**Solves cellular electrophysiology: ion channels, gating variables, ionic currents.**

| Document | Scope |
|----------|-------|
| **[IONIC_MODEL_ARCHITECTURE.md](./ionicModels/IONIC_MODEL_ARCHITECTURE.md)** | Base class, implementations, ODE solvers, tissue variants, GPU acceleration |
| **[README.md](./ionicModels/README.md)** | Available models, registry, single-cell usage |

**Key Models:**
- `ionicModel/` — Abstract base class
- `AlievPanfilov/` — 2-state minimal model
- `BuenoOrovio/` — 4-state reduced model (default, fast)
- `ORd/` — 40+ state comprehensive model (realistic, expensive)
- `TenTusscher/`, `Stewart/`, `Grandi/`, ... — Tissue-specific detailed models
- `*GPU` variants — CUDA implementations for large-scale simulations

**Single-Cell Workflow:**
```bash
singleCellSolver → ionicModel::advance(dt, Vm, gates, Cai)
                → Outputs trace to postProcessing/
                → Compare vs. literature APD, restitution, etc.
```

---

### 3. Verification Models (verificationModels/)

**Validates correctness via manufactured solutions and analytical benchmarks.**

| Document | Scope |
|----------|-------|
| **[VERIFICATION_MODELS_ARCHITECTURE.md](./verificationModels/VERIFICATION_MODELS_ARCHITECTURE.md)** | MMS theory, convergence studies, regression testing |
| **[README.md](./verificationModels/README.md)** | Available tests, how to run |

**Test Categories:**
- `monodomainVerification/` — Monodomain PDE correctness (spatial/temporal order)
- `bidomainVerification/` — Bidomain system (dual potential)
- `ecgVerification/` — Body-surface potential (Laplace equation)

**Typical Usage:**
```bash
# Run on mesh sequence h0, h1, h2, h3
# Compute L2 error for each mesh
# Measure convergence rate p ≈ log(err1/err0) / log(h1/h0)
# Expected: p ≈ 2 (spatial), p ≈ 1–4 (temporal, depending on scheme)
```

---

### 4. Active Tension Models (activeTensionModels/)

**Converts electrical activation → mechanical stress for electro-mechanical coupling.**

| Document | Scope |
|----------|-------|
| **[ACTIVE_TENSION_MODELS_ARCHITECTURE.md](./activeTensionModels/ACTIVE_TENSION_MODELS_ARCHITECTURE.md)** | Phenomenological models, cellular models, parameter fitting |
| **[README.md](./activeTensionModels/README.md)** | Models, configuration, calibration |

**Implementations:**
- `activeTensionModel/` — Abstract base
- `GoktepeKuhl/` — Phenomenological, instantaneous (fast)
- `NashPanfilov/` — Phenomenological, with temporal dynamics (realistic)

**Integration:**
```
Ionic (Vm, Cai) → Active Tension Model (σ_active) → Solid Mechanics (u)
                  ↑___________________________|
                  (feedback: λ_fiber)
```

---

### 5. Coupling Models (couplingModels/)

**Infrastructure for inter-domain coupling (electrical-mechanical, tissue-tissue).**

| Document | Scope |
|----------|-------|
| **[README.md](./couplingModels/README.md)** | Coupler contracts, examples |

**Examples:**
- `electroDomain/` — Domain coupling (Purkinje ↔ myocardium)
- `mechanicalElectroDomain/` — Electro-mechanical feedback
- `fluidSolid/` — Fluid-solid interaction (cardiac hemodynamics, future)

---

### 6. Generic Writer (genericWriter/)

**I/O and post-processing utilities for field output.**

| Document | Scope |
|----------|-------|
| **[README.md](./genericWriter/README.md)** | Output formats, trace generation |

**Features:**
- Volumetric field writing (Vm, Iion, gates, etc.)
- Single-cell trace output (for singleCellSolver, single-point sampling)
- Coordinate-based value queries
- `ionicModelIO` — Formatted output for single-cell data

---

## Quick Navigation: "I want to..."

### "...understand the overall flow"
1. Read **[ARCHITECTURE.md](./ARCHITECTURE.md)** (src-level overview)
2. Read **[electroModels/ELECTROMODEL_ORCHESTRATION.md](./electroModels/core/ELECTROMODEL_ORCHESTRATION.md)**
3. Check **[electroModels/electroCouplers/README.md](./electroModels/electroCouplers/README.md)** for domain coupling

### "...run a monodomain simulation"
1. Read **[electroModels/README.md](./electroModels/README.md)** (configuration overview)
2. Look at **[tutorials/](../tutorials/)** for examples
3. Modify `electroProperties` with `myocardiumSolver monodomainSolver;`
4. Choose `ionicModel` from **[ionicModels/README.md](./ionicModels/README.md)**

### "...validate the solver"
1. Read **[verificationModels/VERIFICATION_MODELS_ARCHITECTURE.md](./verificationModels/VERIFICATION_MODELS_ARCHITECTURE.md)**
2. Run `cd tutorials/manufacturedSolutions/monodomainPseudoECG/` and follow instructions
3. Measure convergence order (should match theory)

### "...add electro-mechanical coupling"
1. Read **[activeTensionModels/ACTIVE_TENSION_MODELS_ARCHITECTURE.md](./activeTensionModels/ACTIVE_TENSION_MODELS_ARCHITECTURE.md)**
2. Set `physicsModel electroMechanicalModel;` in `constant/physicsProperties`
3. Choose `activeTensionModel` and calibrate parameters

### "...implement a new ionic model"
1. Study existing model (e.g., `BuenoOrovio/`) structure
2. Inherit `ionicModel` and implement virtual methods
3. Read **[ionicModels/IONIC_MODEL_ARCHITECTURE.md](./ionicModels/IONIC_MODEL_ARCHITECTURE.md)** Section 12 (Extension Points)
4. Register via `addToRunTimeSelectionTable(ionicModel, MyModel, dictionary)`
5. Test with `singleCellSolver` and compare vs. literature

### "...debug bidirectional Purkinje-myocardium coupling"
1. Read **[electroModels/electroCouplers/README.md](./electroModels/electroCouplers/README.md)** for the PVJ exchange model
2. Read **[electroModels/core/ELECTROMODEL_ORCHESTRATION.md](./electroModels/core/ELECTROMODEL_ORCHESTRATION.md)** for staggered vs PIMPLE scheme behavior
3. Try `electrophysicsAdvanceScheme pimpleStaggeredElectrophysicsAdvanceScheme;`
4. Set `solutionAlgorithm implicit;` to enable PIMPLE iteration
5. Monitor convergence via implicit iterations reported in log

### "...optimize for large meshes"
1. Use GPU ionic model: `ionicModel BuenoOrovioGPU;` or `ORdGPU;`
2. Read **[ionicModels/IONIC_MODEL_ARCHITECTURE.md](./ionicModels/IONIC_MODEL_ARCHITECTURE.md)** Section 7 (GPU Acceleration)
3. Ensure CUDA-capable GPU available and CUDA SDK installed

---

## File Tree with Documentation Pointers

```
src/
├── ARCHITECTURE.md                           # ← You are here (master overview)
├── SRC_DOCUMENTATION_GUIDE.md               # ← This file
│
├── electroModels/
│   ├── ARCHITECTURE.md                      # Domain/solver hierarchy
│   ├── README.md                            # Config overview
│   ├── core/
│   │   ├── README.md                       # Core subsystem overview
│   │   ├── ELECTROMODEL_ORCHESTRATION.md   # ← Core orchestration logic
│   │   ├── electroModel/                    # Entry point
│   │   └── schemes/
│   │       ├── staggered/                   # Weakly coupled (explicit)
│   │       └── pimpleStaggered/             # Strongly coupled (PIMPLE)
│   ├── electroDomains/
│   │   ├── README.md                        # Domain contracts
│   │   ├── myocardiumDomain/                # 3D tissue (primary)
│   │   ├── conductionSystemDomain/          # 1D Purkinje (upstream)
│   │   └── ecgDomain/                       # ECG (downstream)
│   ├── myocardiumModels/
│   │   └── README.md                        # Tissue solver overview
│   ├── conductionSystemModels/
│   │   └── README.md                        # Graph solver overview
│   ├── ecgModels/
│   │   └── README.md                        # ECG solver overview
│   └── electroCouplers/
│       └── README.md                        # Purkinje-myocardium coupling
│
├── ionicModels/
│   ├── IONIC_MODEL_ARCHITECTURE.md          # ← Full ionic model guide
│   ├── README.md                            # Model list, usage
│   ├── ionicModel/                          # Base class
│   ├── BuenoOrovio/                         # 4-state reduced (default)
│   ├── ORd/                                 # 40+ state comprehensive
│   ├── TenTusscher/                         # Detailed ventricular
│   ├── Stewart/                             # Detailed atrial
│   ├── *GPU/                                # CUDA variants
│   └── ...14 more models
│
├── verificationModels/
│   ├── VERIFICATION_MODELS_ARCHITECTURE.md  # ← MMS & convergence guide
│   ├── README.md                            # Available tests
│   ├── monodomainVerification/              # Monodomain MMS
│   ├── bidomainVerification/                # Bidomain MMS
│   ├── ecgVerification/                     # Laplace equation MMS
│   └── ...
│
├── activeTensionModels/
│   ├── ACTIVE_TENSION_MODELS_ARCHITECTURE.md # ← Electro-mechanical guide
│   ├── README.md                            # Models, calibration
│   ├── activeTensionModel/                  # Base class
│   ├── GoktepeKuhl/                         # Phenomenological (fast)
│   └── NashPanfilov/                        # Phenomenological (kinetic)
│
├── couplingModels/
│   ├── README.md                            # Coupler contracts
│   └── common/                              # Shared coupling utilities
│
├── genericWriter/
│   ├── README.md                            # I/O utilities
│   └── ...
│
└── ARCHITECTURE.md                          # Overall src/ structure
```

---

## Key Concepts to Understand

### 1. Runtime Selection Table (RTST)

Models registered at compile-time; selected at runtime via dictionary:

```cpp
// In MyModel.C:
addToRunTimeSelectionTable(ionicModel, BuenoOrovio, dictionary);

// In electroProperties at runtime:
ionicModel  BuenoOrovio;
```

**Benefits:**
- Single binary for all models
- No recompilation to switch models
- Easy to add new models

### 2. Operator Splitting

Electrophysiology simulation splits into two phases:

```
Ionic ODE:      dVm/dt = Iion/Cm + ...
Electrical PDE: ∂Vm/∂t − ∇·(σ∇Vm) = Iion
```

**Advantage:** Use different time steps, solvers, schemes for each without coupling instability.

### 3. Domain Abstraction

Three physical domains connected via **narrow interfaces:**

```
Myocardium    ← provides Vm to couplers
      ↑ ↓    ← receives coupling current
  Coupler     ← couples voltage & current
      ↑ ↓
  Purkinje    ← provides Vp to coupler
```

Each domain is **independent**; couplers are **stateless** (no history).

### 4. Advance Schemes

Different strategies for **when** and **how many times** to call each domain:

- **Explicit staggered:** One pass, fast, weak coupling
- **Implicit pimple:** Multiple passes, convergence, strong coupling

Choose based on coupling strength (Rpvj) and fidelity requirement.

---

## Summary: What Each Layer Does

| Layer | What | Why |
|-------|------|-----|
| **electroModel (core)** | Orchestrates domains & time-stepping | Cleanly separates physics from numerics |
| **Domains** (myocardium, Purkinje, ECG) | Solves PDEs on specific mesh | Each domain independent, easy to extend |
| **Advance Schemes** | Chooses coupling strategy | Fast explicit or stable implicit, user selectable |
| **Ionic Models** | Cellular electrophysiology | Modular, swappable, validated, GPU-capable |
| **Verification Models** | Validates correctness | Finds bugs early via convergence tests |
| **Active Tension** | Electrical → mechanical | Enables electro-mechanical simulation |
| **Couplers** | Inter-domain data exchange | Minimal, focused interfaces |

Together they enable **multi-scale, multi-physics cardiac simulation** from ion channels to whole-heart mechanics.

---

## References

- **Electromodel architecture:** See **[ELECTROMODEL_ORCHESTRATION.md](./electroModels/core/ELECTROMODEL_ORCHESTRATION.md)**
- **Purkinje coupling:** See **[electroCouplers/README.md](./electroModels/electroCouplers/README.md)**
- **Ionic models:** See **[IONIC_MODEL_ARCHITECTURE.md](./ionicModels/IONIC_MODEL_ARCHITECTURE.md)**
- **Verification:** See **[VERIFICATION_MODELS_ARCHITECTURE.md](./verificationModels/VERIFICATION_MODELS_ARCHITECTURE.md)**
- **Mechanics:** See **[ACTIVE_TENSION_MODELS_ARCHITECTURE.md](./activeTensionModels/ACTIVE_TENSION_MODELS_ARCHITECTURE.md)**
