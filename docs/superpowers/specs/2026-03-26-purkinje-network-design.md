# Spec: 1D Purkinje Network Integration

**Date:** 2026-03-26
**Branch:** cardiacFoam+solids4foamv2

---

## Overview

Integrate a graph-based 1D Purkinje Network (PN) solver into the existing 3D biventricular
monodomain framework using the Vergara-Quarteroni operator-splitting coupling approach.
Terminal nodes (PVJs) are pre-aligned to 3D myocardial cell centroids — no spatial search
is required.

---

## Class Structure

Follows the `ecgModel` companion pattern exactly: a separate class hierarchy owned by
`monoDomainElectro`, not a new `electroModel` subclass.

```
purkinjeModel                (abstract base, RTST, analogous to ecgModel)
  └── purkinjeNetworkModel   (concrete, registered via addToRunTimeSelectionTable)
```

`monoDomainElectro` gains one new member:

```cpp
autoPtr<purkinjeModel> purkinjeModelPtr_;   // null if purkinjeNetwork subdict absent
```

### `purkinjeModel` interface

```cpp
class purkinjeModel
{
public:
    TypeName("purkinjeModel");
    declareRunTimeSelectionTable(...);

    purkinjeModel(const volScalarField& Vm, const dictionary& dict);

    static autoPtr<purkinjeModel> New(const volScalarField& Vm, const dictionary& dict);

    virtual ~purkinjeModel() = default;

    virtual void evolve
    (
        scalar t0,
        scalar dt,
        const volScalarField& Vm,
        volScalarField& externalStimulusCurrent
    ) = 0;

    virtual void end() {}
};
```

`evolve()` **adds** the PVJ coupling current into `externalStimulusCurrent` in-place.
No new field and no change to the `solve()` line in `monoDomainElectro` are required.

---

## Dictionary Format

```cpp
// constant/electroProperties
electroModel  monoDomainElectro;

// ... existing monodomain entries ...

purkinjeNetwork
{
    purkinjeModel    purkinjeNetworkModel;

    R_pvj            500.0;    // [Ω·cm²] global PVJ coupling resistance

    // Graph topology: ( nodeA  nodeB  length  conductance )
    // Node 0 is always the root (His bundle / AV node entry)
    edges
    (
        ( 0  1  0.01  0.003 )
        ( 1  2  0.01  0.003 )
        ( 1  3  0.01  0.003 )
    );

    // Terminal (PVJ) nodes and matching 3D cell centroids (parallel lists)
    pvjNodes   ( 2  3 );
    pvjCellIDs ( 42 107 );

    rootStimulus
    {
        startTime    0.0;
        duration     0.002;
        intensity    80000.0;
    }

    purkinjeNetworkModelCoeffs
    {
        ionicModel    Stewart;
    }
}
```

`pvjNodes[i]` and `pvjCellIDs[i]` are parallel — index `i` maps Purkinje terminal node to
3D cell index. All non-root leaf nodes must appear in `pvjNodes`.
`purkinjeNetworkModelCoeffs` follows the `<type>Coeffs` convention of `ecgModel::New`.

---

## Internal State — `purkinjeNetworkModel`

```cpp
// Graph topology (built at construction from 'edges')
label                  nNodes_;
List<List<label>>      neighbours_;      // adjacency list
scalarField            edgeLength_;      // L per edge
scalarField            edgeConductance_; // σ per edge

// PVJ mapping
labelList              pvjNodes_;
labelList              pvjCellIDs_;

// Root stimulus
scalar                 rootStartTime_;
scalar                 rootDuration_;
scalar                 rootIntensity_;

// 1D state — one entry per Purkinje node
scalarField            Vm1D_;
List<scalarField>      ionicStates_;     // [nNodes][NUM_STATES]
List<scalarField>      ionicAlgebraic_;  // [nNodes][NUM_ALGEBRAIC]

// Coupling
scalar                 R_pvj_;

// Ionic model type name (read from dict; drives which CellML functions are called)
word                   ionicModelType_;
```

`ionicModel::solveODE` operates on `volScalarField` arrays and cannot be called per-node
directly. For the first implementation, `ionicModelType_` selects the CellML
`computeVariables` function at construction via a `switch`/factory lookup inside
`purkinjeNetworkModel`. Only `Stewart` is supported in v1; the hook is in place for
extension. The per-node ODE loop iterates over all nodes using the resolved function pointer
— no `volScalarField` wrapper is created.

The explicit 1D graph Laplacian at node `i`:

```
dVm1D[i]/dt = (1/χCm) * Σ_j [ σ_ij*(Vm1D[j]-Vm1D[i]) / L_ij² ] - Iion[i]
```

---

## Call Sequence

### `monoDomainElectro::evolveExplicit()` (modified)

```cpp
// 1) External stimulus (unchanged)
updateExternalStimulusCurrent(externalStimulusCurrent_, externalStimulus_, t0);

// 1b) Advance 1D Purkinje network; adds I_coupling into externalStimulusCurrent_
if (purkinjeModelPtr_)
{
    purkinjeModelPtr_->evolve(t0, dt, Vm_, externalStimulusCurrent_);
}

// 2) Advance 3D ionic model (unchanged)
ionicModelPtr_->solveODE(t0, dt, Vm_, Iion_);

// 3) Solve 3D Vm PDE (unchanged)
solve(chi_*Cm_*fvm::ddt(Vm_) == fvc::laplacian(conductivity_, Vm_)
                               - chi_*Cm_*Iion_
                               + externalStimulusCurrent_);
```

`evolveImplicit()` receives the identical `purkinjeModelPtr_->evolve(...)` call at the same
position.

### `purkinjeNetworkModel::evolve()` steps

1. Apply `rootIntensity_` to `Vm1D_[0]` as a point-source current if `t0` is within the
   root stimulus window.
2. For each node `i`: call CellML `computeVariables` → `Iion1D[i]`.
3. For each node `i`: explicit Euler update using graph Laplacian minus `Iion1D[i]`.
4. For each PVJ `k`:
   ```
   I_coupling = (Vm1D_[pvjNodes_[k]] - Vm[pvjCellIDs_[k]]) / R_pvj_
   externalStimulusCurrent[pvjCellIDs_[k]] += I_coupling
   ```

---

## Files

### New
```
src/electroModels/purkinjeModel/purkinjeModel.H
src/electroModels/purkinjeModel/purkinjeModel.C
src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.H
src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.C
```

### Modified
```
src/electroModels/monoDomainElectro/monoDomainElectro.H   add autoPtr<purkinjeModel>
src/electroModels/monoDomainElectro/monoDomainElectro.C   construct + call in evolve
src/electroModels/Make/files                               add 4 new .C entries
```

No new library dependencies — `ionicModel` and `volFields` are already linked.

---

## Out of Scope (first implementation)

- Interior nodes per branch (each edge = single segment only)
- VTK polydata network input → see `docs/plans/purkinje-explicit-coupling-future.md`
- Explicit shared-current coupling → see `docs/plans/purkinje-explicit-coupling-future.md`
- Per-PVJ heterogeneous resistance
- MPI / parallel decomposition of the 1D network
- `evolveImplicit()` subcycling
