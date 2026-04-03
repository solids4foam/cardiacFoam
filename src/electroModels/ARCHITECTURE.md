# electroModels Architecture Guide

## Overview

The `electroModels` library implements a **multi-domain electrophysiology solver** for cardiac simulations, built on OpenFOAM. The architecture prioritizes **domain separation**, **loose coupling**, and **compositional flexibility**.

---

## Domain Separation Architecture

The system decomposes cardiac electrophysiology into independent domains that evolve together:

```
┌─────────────────────────────────────────────────────────────┐
│                    Myocardium Domain (3D)                    │
│            Reaction-Diffusion Electrophysiology              │
│  - Ionic ODE model (per cell)                                │
│  - Spatial diffusion (Laplacian)                             │
│  - State: Vm (voltage), Iion (ionic current)                │
└─────────────────────────────────────────────────────────────┘
                             ▲
                             │ State Flow (one-way)
                             │
    ┌────────────────────────┴────────────────────────┐
    │                                                 │
┌───────────────┐                             ┌──────────────┐
│  Upstream     │                             │ Downstream   │
│   Domains     │                             │    Domains   │
│               │                             │              │
│ Purkinje 1D   │◄──Coupling Models──►       │  ECG Domain  │
│ Conduction    │   (PVJ coupling)            │              │
│ Network       │                             │ Reads:       │
│               │                             │ - Vm         │
│               │                             │ - σ tensor   │
│               │                             │              │
│ State: Vm     │                             │ No feedback  │
│ per node      │                             │ to myocardium│
└───────────────┘                             └──────────────┘
```

### Domain Categories

1. **Primary Domain**: Myocardium (3D tissue)
   - Solves reaction-diffusion PDE: `∂Vm/∂t = (σ·∇²Vm)/χ + Iion/χ`
   - Owns Vm and activationTime state
   - Exposes state via `electroStateProvider` interface

2. **Upstream Domains**: Optional, advance before myocardium
   - Example: Purkinje network (1D conduction)
   - Can couple bidirectionally with myocardium via coupling models

3. **Downstream Domains**: Optional, advance after myocardium
   - Example: ECG (electro-cardiogram) computation
   - Read-only access to myocardium state

---

## The `electroStateProvider` Interface

**Location**: `core/electroModel/electroStateProvider.H`

### Purpose

`electroStateProvider` is a **read-only state interface** that breaks cyclic dependencies between domains. It enables:
- Downstream domains (ECG) to query myocardium state without providing feedback
- Dependency inversion: domains depend on abstract interface, not concrete solvers
- Graceful handling of optional state channels

### Interface Definition

```cpp
class electroStateProvider
{
public:
    virtual ~electroStateProvider() = default;

    //- Required: Return the computational mesh
    virtual const fvMesh& mesh() const = 0;

    //- Optional: Return transmembrane voltage field (or nullptr)
    virtual const volScalarField* VmPtr() const
    {
        return nullptr;  // Default: not exposed
    }

    //- Optional: Return extracellular potential (bidomain only)
    virtual const volScalarField* phiEPtr() const
    {
        return nullptr;  // Default: not exposed
    }

    //- Optional: Return electrical conductivity tensor
    virtual const volTensorField* conductivityPtr() const
    {
        return nullptr;  // Default: not exposed
    }
};
```

---

## Classes Implementing `electroStateProvider`

### Producers (Expose State)

#### **`electroModel`** (`core/electroModel/electroModel.H`)
- Abstract base class for all EP solvers (monodomain, bidomain, single-cell, eikonal)
- Implements: `Vm()` (pure virtual), optional `VmPtr()`, `phiEPtr()`, `conductivityPtr()`
- Each concrete subclass (e.g., `MonoDomainSolver`) decides which state fields to expose

#### **`MyocardialDomain`** (`domains/MyocardialDomain/MyocardialDomain.H`)
- Abstract base for tissue domains in the domain composition architecture
- Implements: `Vm()` (pure virtual), `sourceField()`, `provider()` for coupling signal access
- Concrete implementation: `ReactionDiffusionMyocardialDomain`

### Consumers (Read State)

#### **`ECGDomain`** (`domains/ECGDomain/ECGDomain.H`)
- Represents ECG electrode measurement domains
- Owns reference to `const electroStateProvider&`
- Queries: `VmPtr()` for voltage, `conductivityPtr()` for anisotropic ECG kernels
- Never writes back; pure read-only consumer

---

## Why `electroStateProvider` Matters

### Problem It Solves

**Without `electroStateProvider`:**
```
ECGModel ──→ MonoDomainSolver ─── bidirectional coupling ──→ [Complex dependencies]
          ──→ BiDomainSolver  ─── hard to extend
          ──→ SingleCellSolver ── can't add new solvers
```

**With `electroStateProvider`:**
```
ECGModel ───→ [abstract] electroStateProvider ───→ MonoDomainSolver
                                               ───→ BiDomainSolver
                                               ───→ Any future solver
                                               (all implement interface)
```

### Benefits

1. **Breaks Cyclic Dependencies**
   - ECG reads myocardium state but doesn't write back
   - Enables independent timestep scheduling
   - Prevents deadlocks in coupled solvers

2. **Extensibility**
   - Add new solver implementations: inherit `electroStateProvider`, automatically compatible with ECG
   - Add new ECG algorithms: query the state provider, works with all solvers

3. **Dependency Inversion**
   - High-level domains depend on abstract interfaces, not concrete implementations

---

## Summary

`electroStateProvider` is the **lynchpin of domain separation** in cardiacFoam:

- **Minimal contract**: Just expose state; don't enforce implementation details
- **One-way flow**: Downstream domains read, never write
- **Composable**: New domains and solvers integrate without modifying existing code
