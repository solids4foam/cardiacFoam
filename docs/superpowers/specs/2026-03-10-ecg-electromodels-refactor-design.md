# ECG ElectroModels Refactor — Design Spec

**Date:** 2026-03-10
**Branch:** `cardiacFoam+solids4foamv2`

---

## Problem

`greensFunctionECGElectro` and `pdeECGElectro` both extend `monoDomainElectro` directly with no shared base, resulting in:

1. **4 copies** of the same conductivity tensor initialisation pattern (`initialiseConductivity`, `initialiseGi` ×2, `initialiseGe`)
2. **Duplicate members** (`Gi_`, `sigmaT_`) declared independently in each subclass
3. **No postProcess mode** — ECG always requires a live monodomain solve; cannot run on existing Vm output
4. **`greensFunctionECGElectro` mixes two distinct ECG approaches** (Green's function BEM and Gima-Rudy dipole) in one class
5. **Dict structure is flat** — monodomain and ECG settings interleaved, no separation

---

## Class Hierarchy

```
electroModel  (abstract, unchanged)
├── monoDomainElectro  (pure Vm solver, unchanged except two additions)
└── ecgElectroBase  (NEW abstract — ECG mode switch + shared state)
    ├── pseudoECGElectro    (Gima-Rudy dipole)
    ├── BEMECGElectro       (Green's function / current source)
    └── pdeECGElectro       (heart Poisson + torso Laplace)
```

`greensFunctionECGElectro` is **removed** — replaced by `pseudoECGElectro` and `BEMECGElectro`.

---

## Layer Responsibilities

| Class | Owns | Does |
|---|---|---|
| `monoDomainElectro` | `conductivity_`, `Vm_`, ionic model | Solves monodomain PDE; adds `readVm()` and `initialiseConductivityTensor()` |
| `ecgElectroBase` | `sigmaI_`, `Gi_`, `sigmaT_`, `postProcess_` | Mode-switches `evolve()`; initialises shared ECG state |
| `pseudoECGElectro` | `gradVm_` | Implements `computeECG()` via Gima-Rudy dipole |
| `BEMECGElectro` | `Is_`, electrodes, torso surface | Implements `computeECG()` via Green's function integral |
| `pdeECGElectro` | `sigmaE_`, `Ge_`, `phiE_`, torso mesh/field | Implements `computeECG()` via Poisson + Laplace |

---

## `monoDomainElectro` Additions

### `initialiseConductivityTensor` (protected static)

Replaces 4 duplicate copies of the disk-or-dict tensor initialisation pattern:

```cpp
static tmp<volTensorField> initialiseConductivityTensor(
    const word& fieldName,
    const dictionary& sourceDict,
    const fvMesh& mesh
);
// 1. Try IOobject READ_IF_PRESENT from disk
// 2. Fallback: read dimensionedSymmTensor from sourceDict, convert via & tensor(I)
```

Used by `monoDomainElectro` for `conductivity_` only. `Gi_` and `Ge_` no longer use this pattern.

### `readVm()` (protected)

```cpp
void readVm();  // calls Vm_.read() — allows ecgElectroBase to update Vm from disk
```

---

## `ecgElectroBase` Interface

### Members

```cpp
protected:
    dimensionedScalar sigmaI_;   // intracellular conductivity scalar (time-invariant)
    volTensorField    Gi_;       // = sigmaI_ * conductivity(), initialised in ctor
    dimensionedScalar sigmaT_;   // isotropic torso / infinite-medium conductivity (constant)
    Switch            postProcess_;
```

`Gi_` and `sigmaT_` are **time-invariant**: initialised once in the constructor, never re-read.

### Virtual interface

```cpp
protected:
    virtual void computeECG() = 0;   // each subclass implements its ECG physics
```

### `evolve()` — not overridable by subclasses

```cpp
void ecgElectroBase::evolve()
{
    if (!postProcess_)
        monoDomainElectro::evolve();   // runs Vm PDE + ionic model
    else
        readVm();                       // reads Vm from current time directory on disk

    computeECG();
}
```

### `read()`

Re-reads `ECGCoeffs` subdict. `sigmaI_` and `sigmaT_` are constants so no field update is needed.

---

## Subclass Physics

### `pseudoECGElectro`

Gima-Rudy dipole approximation:

```
φ_pseudo(x) = (1/sigmaT_) * ∫ (Gi_ & ∇Vm) · (x − x') / |x − x'|³ dV
```

- Members: `gradVm_` (cached, reused each timestep)
- Output: `postProcessing/pseudoECG.dat` (gated on `outputTime()`)

### `BEMECGElectro`

Green's function / boundary element method:

```
Is = −div(Gi_ & ∇Vm)
φ(x) = (1 / 4π·sigmaT_) * ∫ Is / |x − x'| dV
```

- Members: `Is_`, point electrodes, optional torso triSurface
- Output: `postProcessing/ECG.dat` + optional `postProcessing/ECG/phiT_<time>.vtk`
- Parallel-safe via `Foam::reduce(..., sumOp<scalar>())`

### `pdeECGElectro`

Heart Poisson + torso Laplace:

```
div((Gi_ + Ge_) ∇phiE) = −div(Gi_ ∇Vm)   [heart region]
div(sigmaT_ ∇phiT)     = 0                 [torso region]
```

- Members: `sigmaE_` (scalar), `Ge_ = sigmaE_ * conductivity()`, `phiE_`, torso mesh + `phiTPtr_`
- `GiPlusGe` computed inline as `(sigmaI_ + sigmaE_) * conductivity()`
- Output: `phiT.write()` at `outputTime()`

---

## Physical Relationships

```
Gi_ = sigmaI_ * conductivity()   // same anisotropy as monodomain, scaled
Ge_ = sigmaE_ * conductivity()   // same anisotropy as monodomain, scaled
```

`conductivity_` (from `monoDomainElectro`) encodes the full anisotropy tensor. `sigmaI_` and `sigmaE_` are scalar constants read from `ECGCoeffs`.

---

## Dict Structure

```
electroModel  BEMECGElectro;   // or pseudoECGElectro / pdeECGElectro

BEMECGElectroCoeffs
{
    postProcess   false;   // true = skip monodomain solve, read Vm from disk

    monoDomainElectroCoeffs          // ignored when postProcess true
    {
        chi               [0 -1 0 0 0 0 0]    140000;
        Cm                [-1 -4 4 0 0 2 0]   0.01;
        conductivity      [-1 -3 3 0 0 2 0]   (0.17 0 0  0.019 0  0.019);
        ionicModel        BuenoOrovio;
        solutionAlgorithm explicit;
    }

    ECGCoeffs
    {
        sigmaI  [-1 -3 3 0 0 2 0]  0.17;
        sigmaT  [-1 -3 3 0 0 2 0]  0.24725;
        // pdeECGElectro also: sigmaE

        // BEMECGElectro / pseudoECGElectro:
        electrodes
        {
            V1  (-26.1184  -283.543  -70.0558);
        }
        // BEMECGElectro optional:
        // torsoSurface  "constant/triSurface/torso.stl";

        // pdeECGElectro:
        // torsoRegion   torso;
        // heartPatch    epicardium;
        // torsoPatch    endocardium;
    }
}
```

---

## Files Changed

| Action | File |
|---|---|
| Modify | `src/electroModels/monoDomainElectro/monoDomainElectro.H` — add `readVm()`, `initialiseConductivityTensor()` |
| Modify | `src/electroModels/monoDomainElectro/monoDomainElectro.C` — extract utility, add `readVm()` |
| New | `src/electroModels/ecgElectroBase/ecgElectroBase.H` |
| New | `src/electroModels/ecgElectroBase/ecgElectroBase.C` |
| New | `src/electroModels/pseudoECGElectro/pseudoECGElectro.H` |
| New | `src/electroModels/pseudoECGElectro/pseudoECGElectro.C` |
| Rename/refactor | `src/electroModels/greensFunctionECGElectro/` → `BEMECGElectro/` |
| Modify | `src/electroModels/pdeECGElectro/pdeECGElectro.H/.C` — remove duplicate members, extend `ecgElectroBase` |
| Modify | `src/electroModels/Make/files` — add new classes, remove `greensFunctionECGElectro` |
| Modify | `tutorials/ECG/constant/electroProperties` — update dict structure |

---

## Future: Coupled ECG

`coupledECGElectro` will extend `ecgElectroBase` and implement `computeECG()` with bidirectional coupling (extracellular potential feeds back into the monodomain). The base class `evolve()` mode switch accommodates this naturally since `computeECG()` can modify fields used by the next monodomain solve.

---

## Non-Goals

- Anisotropic torso conductivity (future)
- 12-lead ECG derivation from electrode signals (future)
- Multiple monodomain solver types (forward-compatible via `monoDomainElectroCoeffs` subdict naming)
