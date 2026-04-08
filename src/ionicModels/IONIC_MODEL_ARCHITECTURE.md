# Ionic Models — Cellular Electrophysiology Architecture

## 1. Overview

The ionic model layer encapsulates **cellular-scale electrophysiology** — the mathematical models of ion channel kinetics, membrane potential dynamics, and ionic current calculations. Each model represents a specific cell type with specific ion channels and regulatory mechanisms.

This directory provides:
- **ionicModel base class** — abstract interface for all ionic models
- **Concrete implementations** — BuenoOrovio, ORd, TenTusscher, etc. (15+ models)
- **Runtime selection** — switch between models without recompilation
- **Tissue-type variants** — epicardial, endocardial, midmyocardial specializations
- **GPU acceleration** — optional CUDA implementations for high-performance runs

---

## 2. ionicModel Base Class Interface

### Mandatory Virtual Methods

```cpp
class ionicModel
{
    // Advance ionic state by dt using specified ODE solver
    virtual void advance
    (
        scalar dt,
        scalarField& Vm,        // Transmembrane voltage
        scalarField& gates,     // All gating variables (m, h, j, d, f, r, s, ...)
        scalarField& caCyt,     // Intracellular calcium concentration (if tracked)
        const scalarField& Iion // Ionic current from this step
    ) = 0;

    // Return the ionic current given Vm and current gates
    virtual tmp<scalarField> getIionicCurrent
    (
        const scalarField& Vm,
        const scalarField& gates
    ) const = 0;

    // Return number of gating variables (model-dependent)
    virtual label nGates() const = 0;

    // Return equilibrium states (used for initialization)
    virtual tmp<scalarField> initialGates(const scalarField& Vm) const = 0;
    virtual scalar initialVm() const = 0;
};
```

### Optional Virtual Methods

```cpp
    // Compute stimulus current from protocol (S1-S2, pacing, etc.)
    virtual tmp<scalarField> getStimulusCurrent(scalar t) const;

    // Return field names for I/O and post-processing
    virtual wordList exportStateVariables() const;

    // Cell-type specialization selector
    virtual void setTissueType(const word& tissue) = 0;
```

---

## 3. Concrete Model Implementations

### Representative Models

| Model | Type | Complexity | Use Case |
|-------|------|-----------|----------|
| **AlievPanfilov** | Minimal | 2 states (Vm, h) | Fast prototype, educational |
| **BuenoOrovio** | Reduced | 4 states (Vm, u, v, w) | Standard cardiac research |
| **TenTusscher** | Detailed | 19 states | Ventricular cell kinetics |
| **ORd (O'Hara-Rudy)** | Comprehensive | 40+ states | Most realistic; high cost |
| **Stewart** | Detailed | 20+ states | Atrial-specific |

### Model Structure Example: BuenoOrovio

```cpp
class BuenoOrovio : public ionicModel
{
    // Model parameters (tissue-dependent)
    scalar tau_v1, tau_v2, tau_w1, tau_w2, ...;

    // Equilibrium relationships (gating steady-states)
    scalar v_infty(scalar Vm) const;
    scalar w_infty(scalar Vm) const;

    // Time constants (voltage-dependent)
    scalar tau_v(scalar Vm) const;
    scalar tau_w(scalar Vm) const;

    // Ionic current formulation
    // Iion = (Vm - E_K) * g_K * w + (Vm - E_Ca) * g_Ca * d
    //      + leak currents + stimulus

    scalar nGates() const override { return 3; }  // {u, v, w}
};
```

---

## 4. Tissue Type Specialization

Most models support multiple cell types with different kinetics:

```cpp
ionicModel::setTissueType(const word& tissue);
```

**Available specializations** (model-dependent):

| Type | Location | Characteristics |
|------|----------|-----------------|
| **epicardial** | Heart surface | Faster repolarization, shorter APD |
| **endocardial** | Inner wall | Longer APD, more calcium handling |
| **midmyocardial** | Mid-wall | Intermediate properties |
| **atrial** | Atrial tissue | Different ion channel expression |
| **ventricular** | Ventricular tissue | Standard single cell |

Example: BuenoOrovio specializations change gate time constants `tau_v`, `tau_w`, etc.

---

## 5. ODE Time Integration Strategies

### Explicit Methods (Fast, May Require Small Δt)

```cpp
explicit:
{
    solver RKF45;          // Runge-Kutta-Fehlberg 4(5)
    solver RK4;            // Runge-Kutta 4th order (fixed step)
    solver euler;          // Forward Euler (rarely used; unstable)

    initialODEStep 1e-6;   // First substep for adaptive solver
    maxSteps 1000000000;   // Hard limit to prevent infinite loops
    absTol 1e-6;           // Local error tolerance
    relTol 1e-4;           // Relative error tolerance
}
```

**RKF45 (Default):**
- Adaptive step control: accepts/rejects substeps based on error estimate
- Efficient for stiff and non-stiff problems
- Automatic step size adjustment ensures accuracy

### Implicit Methods (Stable for Large Δt; Higher Cost)

```cpp
implicit:
{
    solver CVODE;          // SUNDIALS ODE solver (requires PETSc/CVODE)
    solver DASPK;          // DAE solver

    // Implicit solvers are externally integrated; not shown here
}
```

**Trade-off:**
- Explicit RKF45 + small Δt: stable, faster per step, many steps needed
- Implicit + large Δt: more cost per step, fewer steps needed, overall faster for stiff problems

---

## 6. Input/Output and State Management

### State Vector Layout

Each cell has a state vector: `[Vm, gate1, gate2, ..., gateN, caCyt, ...]`

```cpp
label stateIndex = 0;  // Vm
label gatesStart = 1;  // First gating variable
label caCytIndex = 1 + nGates();  // Intracellular Ca (if model uses it)
```

### Export Variables

Selected state variables are written to volumetric fields:

```cpp
// In electroProperties:
outputVariables
{
    ionic
    {
        export  (Vm Iion gates);   // Written to disk
        debug   (Vm);              // Printed to terminal each step
    }
}
```

### Single-Cell Traces

For validation/fitting:

```cpp
// In singleCellSolver runs:
singleCellStimulus
{
    stim_start      20;      // [ms] stimulus onset
    stim_duration   1;       // [ms]
    stim_amplitude  0.4;     // [nA/pF]
    stim_period_S1  1000;    // [ms] pacing period
    nstim1          3;       // number of S1 beats
    stim_period_S2  200;     // [ms] S2 timing (0 = no S2)
    nstim2          1;       // number of S2 beats
}
```

Output: `postProcessing/<ionicModel>_<tissue>_<stimulus>.txt` with full state vector time series.

---

## 7. GPU Acceleration

Selected models have CUDA implementations:

```cpp
BuenoOrovioGPU
ORdGPU, ORdGPUCompact, ORdGPUOpt
TNNPGPU
```

**When to use:**
- ✅ Large 3D meshes (>1M cells)
- ✅ Long simulations (many time steps)
- ⚠️ Requires NVIDIA GPU with CUDA capability
- ⚠️ GPU memory constraints for very large models

**Trade-off:**
- Excellent for forward simulations (batch compute)
- Less ideal for inverse problems (gradient computation not yet available)

---

## 8. Stimulus Protocol Architecture

### Volume Stimulus (Monodomain PDE Level)

```cpp
monodomainStimulus
{
    stimulusLocationMin  (0 0 5.5e-3);   // Stimulus box corner 1 [m]
    stimulusLocationMax  (1.5e-3 1.5e-3 7e-3);  // Stimulus box corner 2
    stimulusDuration     [0 0 1 0 0 0 0]  2e-3;  // [s]
    stimulusIntensity    [0 -3 0 0 0 1 0] 50000; // [A/m^3]
    stimulusStartTime    0.0;            // [s]
}
```

**Applied as:** volumetric source term in myocardium domain PDE.

### Single-Cell Stimulus (Ionic Model Level)

```cpp
singleCellStimulus
{
    stim_amplitude  0.4;    // [nA/pF] current per unit membrane capacitance
    // Converts to current density for single cell:
    // Istim = 0.4 * (membrane area) / (cell volume) [A/m^3]
}
```

**Applied as:** single point current in ionic model ODE.

### Pacing Protocols

Both volume and cell-level stimuli support complex protocols:
- **S1 fixed pacing**: Regular beats at interval S1
- **S1-S2 protocol**: N beats at S1, then one beat at S2 (restitution, vulnerability studies)
- **Custom**: User-defined time-varying `Istim(t)`

---

## 9. Initialization and Equilibration

### Resting State

Most models define a resting transmembrane voltage (e.g., -84 mV):

```cpp
scalar Vm_rest = ionicModel->initialVm();  // e.g., -0.084 V
```

At rest, all gates are at equilibrium:

```cpp
scalarField gates_rest = ionicModel->initialGates(Vm_rest);
```

This is the default IC for transient simulations.

### Equilibration Runs

Before main simulation, users often run a "steady-state pacing" phase:
- Apply S1 stimulus repeatedly (e.g., 100 beats at 1 Hz)
- Let ionic state converge to periodic steady-state
- Use final state as IC for transient experiment

This accounts for slow secondary variables (Ca2+ accumulation, channel slow inactivation) that don't reach equilibrium in one beat.

---

## 10. Key Design Decisions

### Why Runtime Selection?

- Users fit models to different datasets (ORd for normal, Stewart for atrial, etc.)
- Hardware constraints (GPU availability) drive model choice
- Code compiled once; model swapped via dictionary

### Why Separate ODE from PDE?

- Decouples cellular kinetics from tissue-scale physics
- Allows independent testing/validation of ionic model
- Enables operator-splitting schemes (e.g., Godunov, Strang)

### Why Tissue Types Over Separate Classes?

- Avoids 15 models × 3 cell types = 45 separate C++ classes
- Kinetics differ slightly; same structure
- Runtime switch via `setTissueType()` is clean and flexible

---

## 11. Common Integration Patterns

### Pattern 1: Electrokinetic Coupling via Monodomain

```
Voltage PDE:  ∂Vm/∂t − ∇·(σ∇Vm) = −Iion + Istim
Ionic ODE:    dgate/dt = f(Vm, gate)
              Iion = h(Vm, gates)
```

At each cell:
1. Receive Vm from PDE (from previous time step or iteration)
2. Advance ionic ODE: `gate(t+Δt) = advance(gate(t), Vm, Δt)`
3. Compute `Iion(t+Δt)` from updated gates
4. Return to PDE as source term

### Pattern 2: Single-Cell Benchmark

```
ODE only:     dVm/dt = −Iion / Cm + Istim / Cm
              dgate/dt = f(Vm, gate)
```

No spatial coupling; all cells independent. Used for:
- Model validation against literature
- Restitution curve measurement
- Action potential duration fitting

---

## 12. Extension Points

### Adding a New Ionic Model

1. Inherit `ionicModel`
2. Implement all pure virtual methods
3. Define state vector layout (Vm index, gate indices, Ca index, etc.)
4. Implement advance logic (call RKF45 or similar ODE integrator)
5. Register: `addToRunTimeSelectionTable(ionicModel, MyModel, dictionary)`
6. Document tissue types, model parameters, reference paper

### Adding Tissue Type Specialization

1. Add parameters for new tissue variant
2. Implement conditional logic in constructors/setters based on tissue name
3. Update `setTissueType()` to recognize new tissue
4. Test equilibrium action potentials match literature

---

## 13. Summary

Ionic models are the **cellular foundation** of cardiac electrophysiology:

| Aspect | Purpose | Example |
|--------|---------|---------|
| **Base class** | Common interface | `ionicModel::advance()` |
| **Concrete models** | Cellular kinetics | BuenoOrovio (4 gates), ORd (40+ states) |
| **ODE solvers** | Time integration | RKF45 (adaptive), CVODE (implicit) |
| **Tissue types** | Cell-type variability | epicardial, endocardial specializations |
| **Stimulus** | Excitation protocols | pacing, S1-S2, volumetric/point |
| **GPU acceleration** | High performance | CUDA variants for big simulations |

The modular design allows researchers to seamlessly combine models with different fidelity levels, tissue representations, and computational backends in a single simulation.

