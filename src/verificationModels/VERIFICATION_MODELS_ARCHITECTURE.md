# Verification Models — Manufactured Solutions and Validation Infrastructure

## 1. Overview

The verification models layer provides **analytical and manufactured solutions** to validate the correctness of the solver infrastructure. These models replace the ionic kinetics with **exact (or analytically known) solutions**, allowing bit-for-bit comparison between computed results and theoretical expectations.

This is critical because:
- Ionic models are complex nonlinear systems (difficult to verify analytically)
- Manufactured solutions have **zero model error** — all error is discretization error
- Convergence studies isolate spatial/temporal order of accuracy
- Regression tests catch regressions before they reach production

---

## 2. Categories of Verification Models

### A. Manufactured (MMS) Models

Create a **smooth, spatially-varying analytical Vm(x,y,z,t)** that satisfies the PDE exactly, then reverse-engineer the required source term.

**Example: ManufacturedFDA**

Vm is manufactured as a separable product:

```
Vm(x,y,z,t) = sin(πx/Lx) × sin(πy/Ly) × sin(πz/Lz) × sin(πt/T)
```

This smooth function:
- Has known spatial derivatives (∇Vm, ∇²Vm)
- Has known time derivative (∂Vm/∂t)
- Satisfies any linear PDE exactly if we add the right source term

For monodomain:
```
∂Vm/∂t − ∇·(σ∇Vm) = S_mfg(x,y,z,t)
```

We compute: `S_mfg = ∂Vm/∂t − σ∇²Vm` analytically, then inject it as the RHS.

**Advantages:**
- ✅ Exact solution everywhere in domain
- ✅ Arbitrary convergence order (choose polynomial degree)
- ✅ Works for both explicit and implicit schemes
- ⚠️ Ignores ionic model physics (no Iion feedback)

### B. Manufactured Stimulus (Non-Ionic) Models

Replace ionic currents with **manufactured stimulus patterns** while keeping PDE structure intact.

**Example: monodomainFDAManufactured**

```
Iion_manufactured(x,y,z,t) = A × sin(ωx × x) × sin(ωy × y) × sin(ωt × t)
```

**Characteristics:**
- ✅ Includes spatial coupling (∇·(σ∇Vm) term is real diffusion)
- ✅ Tests ionic integration + PDE solve together
- ⚠️ Still lacks true ionic model feedback

### C. Null-Ionic Models (Passive Tissue)

Set `Iion ≡ 0` and drive with volumetric stimulus only.

**Use case:**
- Validates diffusion discretization in isolation
- Easiest to achieve high convergence order
- Baseline before testing with real ionic models

---

## 3. Verification Model Implementations

### monodomainVerification

Validates the **monodomain PDE** (single-potential diffusion):

```cpp
class tManufacturedFDAMonodomainVerifier : public ionicModel
{
    // Instead of ionic state variables, use manufactured stimulus
    scalar getStimulusCurrent(scalar t) const override
    {
        return A_ * sin(omega_x_ * x_) * sin(omega_t_ * t_);
    }

    // Advance: just return zero (no ionic dynamics)
    void advance(...) override { /* no-op */ }

    // Provide error norm relative to analytical solution
    scalar L2_error(const scalarField& computed_Vm) const;
};
```

### bidomainVerification

Validates the **bidomain system** (dual potential):

```cpp
∂Vm/∂t − ∇·(σi∇Vm) − ∇·(σi∇φe) = −Iion
0 = ∇·(σi + σe)∇φe + ∇·(σi∇Vm)
```

More complex than monodomain:
- Two unknowns (Vm, φe)
- Algebraic constraint on φe (no time derivative)
- Requires specialized preconditioners for implicit solve

**Verification approach:**
- Manufacture `Vm(x,y,z,t)` (as before)
- Compute required σi and σe tensors to satisfy both equations
- Validate that computed solution converges to manufactured solution

### ecgVerification

Validates the **ECG domain** in isolation:

```cpp
∇·(σbody ∇φ) = 0  (Laplace eq in body)

Boundary: φ = −Vm(x)  (on heart surface, from myocardium)
```

**Verification:**
- Use manufactured Vm on heart boundary
- Compute φ analytically in simple geometries (sphere, ellipsoid)
- Compare computed ECG potentials at electrode sites

---

## 4. Convergence Study Workflow

### Step 1: Run Simulations on Refined Mesh Sequence

```
Mesh 0: h₀ (coarse, e.g., 2 cm elements)
  ↓ (half element size)
Mesh 1: h₁ = h₀/2
  ↓ (half again)
Mesh 2: h₂ = h₁/2
  ↓ (half again)
Mesh 3: h₃ = h₂/2 (fine, e.g., 0.25 cm elements)
```

For each mesh:
- Run with manufactured solution
- Compute L2 error: `||computed_Vm − analytical_Vm||_L2`
- Record error and mesh size h

### Step 2: Compute Convergence Rate

With results (h, error):

```
log(error) ≈ log(C) + p × log(h)
```

Least-squares fit gives **convergence order p**:

```
p = [log(err₁/err₀)] / [log(h₁/h₀)]
```

**Expected results:**

| Scheme | Space Order | Time Order | Combined |
|--------|-------------|-----------|----------|
| FVM + linear basis | 2 | 1 | O(h² + Δt) |
| Explicit RK4 + FVM | 2 | 4 | O(h² + Δt⁴) |
| IMEX scheme + FVM | 2 | 2 | O(h² + Δt²) |

If computed order matches theory → **code is correct** (discretization is as intended).

### Step 3: Regression Testing

Establish a baseline on coarse mesh:

```
baseline_error = 1.234e-3  (for h = 2 cm)
```

In continuous integration:
- Run on same coarse mesh
- If `new_error > 1.5 × baseline_error` → **FAIL** (regression detected)
- If within tolerance → **PASS**

This catches bugs before they spread to users.

---

## 5. Configuration: electroProperties Setup

### Using Monodomain Verification

```cpp
monodomainSolverCoeffs
{
    // Replace actual ionic model with test model
    ionicModel  monodomainFDAManufactured;    // Manufactured stimulus, no ionic dynamics
    tissue      myocyte;              // Tissue type (usually ignored for MMS)

    // Conductivity (choose something reasonable)
    conductivity [-1 -3 3 0 0 2 0] (0.1334 0 0 0.01761 0 0.01761);
    chi   [0 -1 0 0 0 0 0]  140000;
    cm    [-1 -4 4 0 0 2 0]  0.01;

    // No monodomainStimulus block needed; verifier provides Iion(t)
}
```

### Single-Cell Verification

```cpp
// In singleCellSolver:
ionicModel  monodomainFDAManufactured;

singleCellStimulus
{
    // Ignored; verifier provides its own Istim(t)
}
```

Writes single-cell trace to `postProcessing/` and compares against analytical solution point-wise.

---

## 6. Code Integration: How Verifiers Hook In

### Constructor Pattern

```cpp
class tManufacturedFDAMonodomainVerifier : public ionicModel
{
    // Store mesh and time info
    const fvMesh& mesh_;
    const Time& runTime_;

    // Manufactured solution parameters
    scalar A_, omega_x_, omega_y_, omega_z_, omega_t_;
    scalar Lx_, Ly_, Lz_, T_;

public:
    // Constructor reads verification parameters from dictionary
    tManufacturedFDAMonodomainVerifier
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const Time& runTime
    );

    // Provide manufactured Iion at current time
    scalar getStimulusCurrent(scalar t) const override
    {
        return A_ * sin(omega_x_ * x_global_) * sin(omega_t_ * t);
    }

    // After simulation, compute error norms
    void postProcess(const scalarField& computed_Vm)
    {
        scalarField analytical_Vm = computeAnalyticalSolution();
        scalar L2_error = computeL2Norm(computed_Vm − analytical_Vm);
        writeErrorReport(L2_error);
    }
};
```

### Integration Points

1. **Initialization:** Read manufactured parameters from dict
2. **Each ionic step:** Return manufactured stimulus instead of ionic kinetics
3. **Output phase:** Compute and log error norms
4. **Postprocessing:** Generate convergence plots

---

## 7. Key Verification Tests in Repository

### Available Tests

| Test | File | Purpose |
|------|------|---------|
| **Monodomain MMS** | `monodomainVerification/` | Spatial/temporal convergence of mono PDE |
| **Bidomain MMS** | `bidomainVerification/` | Bidomain system correctness |
| **ECG MMS** | `ecgVerification/` | Body-surface potential (Laplace eq) |
| **Electrode Verification** | (implicit) | ECG at specific electrode sites |

### Running a Verification Test

```bash
cd tutorials/manufacturedSolutions/monodomainPseudoECG/
python prepare_simulation_data.py --refine 4  # Create 4 mesh levels
./Allwmake  # Compile with verification models enabled

# Run coarse and fine
cardiacFoam -case case_h0 > log 2>&1
cardiacFoam -case case_h3 > log 2>&1

# Compare errors
python convergence_analysis.py  # Fit convergence rate
```

Output: ASCII convergence plot showing slope (should match theory).

---

## 8. Limitations and Extensions

### Current Limitations

- ✗ Manufactured solutions often **artificially smooth** (no shock-like activation)
- ✗ Cannot test ionic model directly (verifiers disable ionic dynamics)
- ✗ Requires mesh/time refinement studies (computationally expensive for 3D)

### Future Extensions

1. **Hybrid Verification:** Part of domain uses MMS, rest uses real ionic model
   - Validates coupling without global analytical solution

2. **Adjoint-Based Verification:** Run adjoint equation to assess sensitivity
   - Detects instabilities in non-linear terms

3. **Empirical Convergence:** Compare against high-precision reference solution (ORd on ultra-fine mesh)
   - Avoids manufactured solution bias
   - More representative of real problem

---

## 9. Key Design Decisions

### Why Separate Verifiers?

- Modular: each verifier is a small ionicModel subclass
- Avoids bloating production code with verification logic
- Different models test different aspects (mono vs bi vs ECG)

### Why Manufactured vs Empirical?

- **Manufactured:** Exact derivatives, clear error analysis
- **Empirical:** Realistic physics, but requires reference solution generation

Both are complementary; neither alone is sufficient.

### Why Convergence Studies?

- Single-grid error is **meaningless** without knowing the reference solution
- Convergence rate **proves correctness of implementation**
  - If order drops → implementation bug
  - If matches theory → code is (likely) correct

---

## 10. Summary

Verification models are the **correctness guardrails** of cardiacFoam:

| Component | Purpose | Example |
|-----------|---------|---------|
| **Manufactured solution** | Known analytical Vm(x,y,z,t) | sin products, separable |
| **Convergence study** | Measure discretization error | h refinement → error plot |
| **Regression test** | Catch code changes | CI baseline monitoring |
| **Isolate physics** | Test parts independently | Monodomain (no coupling), ECG (passive) |

They ensure that numerical results are **trustworthy** before clinical deployment or publication.
