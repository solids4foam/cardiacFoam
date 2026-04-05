# ECG ElectroModels Refactor — Design Spec

**Date:** 2026-03-10
**Branch:** `cardiacFoam+solids4foamv2`

---

## Problem

`greensFunctionECGElectro` and `pdeECGElectro` both extend `MonoDomainSolver` directly with no shared base, resulting in:

1. **4 copies** of the same conductivity tensor initialisation pattern
2. **Duplicate members** (`Gi_`, `sigmaT_`) declared independently in each subclass
3. **No postProcess mode** — ECG always requires a live monodomain solve
4. **`greensFunctionECGElectro` mixes two distinct ECG approaches** in one class
5. **No runtime-selectability of ECG type** — changing the ECG approach requires changing `electroModel`

---

## Class Hierarchy

```
electroModel  (abstract, unchanged)
├── MonoDomainSolver  (pure Vm solver — unchanged except two additions)
├── ecgElectro         (NEW concrete — owns autoPtr<ecgModel>, mode switch)
└── pdeECGElectro      (extends MonoDomainSolver — stays separate, manages torso mesh)

ecgModel  (NEW abstract, runtime-selectable — same pattern as ionicModel)
├── BEMECGElectro      (Green's function / current source)
└── pseudoECGElectro   (Gima-Rudy dipole)
```

`greensFunctionECGElectro` is **removed** — replaced by `ecgElectro` + `BEMECGElectro`/`pseudoECGElectro`.

`pdeECGElectro` stays as a direct `MonoDomainSolver` subclass because it manages a separate torso mesh and solves additional PDEs — too coupled for the `ecgModel` pattern.

---

## Layer Responsibilities

| Class | Owns | Does |
|---|---|---|
| `MonoDomainSolver` | `conductivity_`, `Vm_`, ionic model | Solves monodomain PDE; adds `readVm()` and `initialiseConductivityTensor()` |
| `ecgElectro` | `postProcess_`, `autoPtr<ecgModel>` | Mode-switches `evolve()`; delegates ECG computation to `ecgModel` |
| `ecgModel` (abstract) | `Gi_`, `sigmaT_` | Defines runtime-selectable ECG compute interface |
| `BEMECGElectro` | `Is_`, electrodes, torso surface | Implements `compute()` via Green's function integral |
| `pseudoECGElectro` | electrodes | Implements `compute()` via Gima-Rudy dipole |
| `pdeECGElectro` | `Gi_`, `Ge_`, `sigmaT_`, `phiE_`, torso mesh | Solves heart Poisson + torso Laplace; not an `ecgModel` |

---

## `MonoDomainSolver` Additions

### `initialiseConductivityTensor` (protected static)

Replaces duplicate tensor initialisation pattern in both monodomain and ECG subclasses:

```cpp
static tmp<volTensorField> initialiseConductivityTensor(
    const word& fieldName,
    const dictionary& sourceDict,
    const fvMesh& mesh
);
// 1. Try IOobject READ_IF_PRESENT from disk
// 2. Fallback: read dimensionedSymmTensor from sourceDict, convert via & tensor(I)
```

Used by `MonoDomainSolver` (for `conductivity_`) and by `ecgModel` subclasses (for `Gi_`). Since it is `static`, `ecgModel` subclasses call it as `MonoDomainSolver::initialiseConductivityTensor(...)` after including the header via `lnInclude/`.

### `readVm()` (protected)

```cpp
void readVm();  // calls Vm_.read() — used by ecgElectro in postProcess mode
```

---

## `ecgModel` Interface

New abstract class, follows the same runtime-selection pattern as `ionicModel`.

```cpp
class ecgModel
{
protected:
    const volScalarField& Vm_;   // reference to monodomain Vm (owned by ecgElectro)
    const fvMesh&         mesh_;

    volTensorField        Gi_;   // intracellular conductivity (time-invariant)
    dimensionedScalar     sigmaT_; // torso / infinite-medium conductivity

public:
    declareRunTimeSelectionTable(...);

    static autoPtr<ecgModel> New(
        const volScalarField& Vm,
        const dictionary& dict      // the ecgModel's own Coeffs subdict
    );

    virtual void compute() = 0;    // called every timestep after Vm is ready
    virtual bool read()   = 0;
};
```

`Gi_` is initialised via `MonoDomainSolver::initialiseConductivityTensor("Gi", dict, mesh)`. `Gi_` and `sigmaT_` are **time-invariant** — set once in the constructor.

---

## `ecgElectro` Interface

Concrete `electroModel`, replaces `ecgElectroBase` as the user-facing electroModel name.

```cpp
class ecgElectro : public MonoDomainSolver
{
    const Switch postProcess_;
    autoPtr<ecgModel> ecgModelPtr_;

public:
    TypeName("ecgElectro");

    virtual bool evolve()
    {
        bool ok = true;
        if (!postProcess_)
            ok = MonoDomainSolver::evolve();
        else
            readVm();

        ecgModelPtr_->compute();
        return ok;
    }
};
```

---

## Subclass Physics

### `BEMECGElectro` (ecgModel subclass)

Green's function / boundary element method:

```
Is = −div(Gi_ & ∇Vm)
φ(x) = (1 / 4π·sigmaT_) * ∫ Is / |x − x'| dV
```

- Members: `Is_`, point electrodes, optional torso triSurface
- Output: `postProcessing/ECG.dat` + optional `postProcessing/ECG/phiT_<time>.vtk`

### `pseudoECGElectro` (ecgModel subclass)

Gima-Rudy dipole approximation:

```
φ_pseudo(x) = (1/sigmaT_) * ∫ (Gi_ & ∇Vm) · (x − x') / |x − x'|³ dV
```

- Members: point electrodes
- Output: `postProcessing/pseudoECG.dat`

### `pdeECGElectro` (MonoDomainSolver subclass — unchanged pattern)

Heart Poisson + torso Laplace — stays as a direct `electroModel` subclass because it needs full PDE infrastructure on a separate torso mesh. Retains its own `Gi_`, `Ge_`, `sigmaT_` for now (Gi/Ge tensor → scalar simplification deferred).

---

## Dict Structure

```
electroModel  ecgElectro;

ecgElectroCoeffs
{
    postProcess   false;   // true = skip monodomain solve, read Vm from disk

    // Monodomain settings (flat — same keys as MonoDomainSolverCoeffs):
    chi               [0 -1 0 0 0 0 0]    140000;
    Cm                [-1 -4 4 0 0 2 0]   0.01;
    conductivity      [-1 -3 3 0 0 2 0]   (0.17 0 0  0.019 0  0.019);
    ionicModel        BuenoOrovio;
    solutionAlgorithm explicit;
    // ... other monodomain keys unchanged ...

    // ECG type (runtime-selectable, same pattern as ionicModel):
    ecgModel      BEMECGElectro;   // or pseudoECGElectro

    BEMECGElectroCoeffs
    {
        Gi      [-1 -3 3 0 0 2 0] (0.17 0 0  0.019 0  0.019);
        sigmaT  [-1 -3 3 0 0 2 0] 0.24725;
        electrodes
        {
            V1  (-26.1184  -283.543  -70.0558);
            // ...
        }
        // optional: torsoSurface "constant/triSurface/torso.stl";
    }

    pseudoECGElectroCoeffs
    {
        Gi      [-1 -3 3 0 0 2 0] (0.17 0 0  0.019 0  0.019);
        sigmaT  [-1 -3 3 0 0 2 0] 0.24725;
        electrodes { ... }
    }
}
```

---

## Files Changed

| Action | File |
|---|---|
| Modify | `src/electroModels/MonoDomainSolver/MonoDomainSolver.H` — add `readVm()`, `initialiseConductivityTensor()` |
| Modify | `src/electroModels/MonoDomainSolver/MonoDomainSolver.C` — implement both |
| **New** | `src/electroModels/ecgModel/ecgModel.H` — abstract ecgModel base + run-time selector |
| **New** | `src/electroModels/ecgModel/ecgModel.C` |
| **New** | `src/electroModels/ecgElectro/ecgElectro.H` — concrete electroModel wrapper |
| **New** | `src/electroModels/ecgElectro/ecgElectro.C` |
| **New** | `src/electroModels/ecgModels/BEMECGElectro/BEMECGElectro.H` |
| **New** | `src/electroModels/ecgModels/BEMECGElectro/BEMECGElectro.C` |
| **New** | `src/electroModels/ecgModels/pseudoECGElectro/pseudoECGElectro.H` |
| **New** | `src/electroModels/ecgModels/pseudoECGElectro/pseudoECGElectro.C` |
| Modify | `src/electroModels/pdeECGElectro/pdeECGElectro.H/.C` — unchanged hierarchy, cleanup only |
| Delete | `src/electroModels/greensFunctionECGElectro/` |
| Modify | `src/electroModels/Make/files` — add new classes, remove `greensFunctionECGElectro` |
| Modify | `tutorials/ECG/constant/electroProperties` — `electroModel ecgElectro`, add `ecgModel BEMECGElectro` |

---

## Future: Coupled ECG

A future `coupledECGElectro` ecgModel subclass will implement bidirectional coupling. The `ecgElectro::evolve()` mode switch already supports this since `ecgModel::compute()` can write fields read by the next monodomain iteration.

---

## Non-Goals (this refactor)

- `Gi`/`Ge` tensor → scalar simplification (deferred)
- Anisotropic torso conductivity
- 12-lead ECG derivation
- `postProcess` time-looping solver integration (flag is wired, looping deferred)
