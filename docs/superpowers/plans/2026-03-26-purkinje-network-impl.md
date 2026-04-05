# Purkinje Network Integration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a graph-based 1D Purkinje Network (PN) solver that couples into `MonoDomainSolver` via operator splitting, injecting PVJ coupling currents into the existing `externalStimulusCurrent_` field.

**Architecture:** A `purkinjeModel` companion class (like `ecgModel`) is owned by `MonoDomainSolver`. The concrete `purkinjeNetworkModel` reads a graph edge list from the dict, advances a 1D explicit-Euler monodomain on the graph each timestep, and writes coupling currents at PVJ cell indices into `externalStimulusCurrent` before the 3D PDE solve.

**Tech Stack:** OpenFOAM C++, wmake, existing `ionicModel` RTST (Stewart_2009 for Purkinje cells), `volScalarField`/`scalarField`.

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `src/electroModels/purkinjeModel/purkinjeModel.H` | Abstract base, RTST declaration |
| Create | `src/electroModels/purkinjeModel/purkinjeModel.C` | RTST definitions, `New` selector |
| Create | `src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.H` | Concrete class: graph, ionic model, state |
| Create | `src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.C` | Constructor + `evolve()` |
| Modify | `src/electroModels/MonoDomainSolver/MonoDomainSolver.H` | Add `autoPtr<purkinjeModel>` |
| Modify | `src/electroModels/MonoDomainSolver/MonoDomainSolver.C` | Construct + call in evolveExplicit/Implicit |
| Modify | `src/electroModels/Make/files` | Add 2 new `.C` source entries |

No changes to `Make/options` — `ionicModels` and `volFields` are already linked.

---

## Task 1: Abstract base `purkinjeModel`

**Files:**
- Create: `src/electroModels/purkinjeModel/purkinjeModel.H`
- Create: `src/electroModels/purkinjeModel/purkinjeModel.C`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1.1: Create `purkinjeModel.H`**

```cpp
// src/electroModels/purkinjeModel/purkinjeModel.H
#ifndef purkinjeModel_H
#define purkinjeModel_H

#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class purkinjeModel
{
protected:
    const volScalarField& Vm_;
    const fvMesh&         mesh_;

public:
    TypeName("purkinjeModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        purkinjeModel,
        dictionary,
        (
            const volScalarField& Vm,
            const dictionary&     dict,
            const scalar          initialDeltaT
        ),
        (Vm, dict, initialDeltaT)
    );

    purkinjeModel
    (
        const volScalarField& Vm,
        const dictionary&     dict
    );

    static autoPtr<purkinjeModel> New
    (
        const volScalarField& Vm,
        const dictionary&     dict,
        const scalar          initialDeltaT
    );

    virtual ~purkinjeModel() = default;

    // Advance the 1D network one timestep and add coupling currents
    // into externalStimulusCurrent at PVJ cell indices.
    virtual void evolve
    (
        scalar                t0,
        scalar                dt,
        const volScalarField& Vm,
        volScalarField&       externalStimulusCurrent
    ) = 0;

    virtual void end() {}
};

} // End namespace Foam

#endif
```

- [ ] **Step 1.2: Create `purkinjeModel.C`**

```cpp
// src/electroModels/purkinjeModel/purkinjeModel.C
#include "purkinjeModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(purkinjeModel, 0);
defineRunTimeSelectionTable(purkinjeModel, dictionary);

purkinjeModel::purkinjeModel
(
    const volScalarField& Vm,
    const dictionary&     dict
)
:
    Vm_(Vm),
    mesh_(Vm.mesh())
{}

autoPtr<purkinjeModel> purkinjeModel::New
(
    const volScalarField& Vm,
    const dictionary&     dict,
    const scalar          initialDeltaT
)
{
    const word modelType(dict.lookup("purkinjeModel"));

    Info<< "Selecting purkinjeModel " << modelType << nl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown purkinjeModel type " << modelType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<purkinjeModel>
    (
        ctorPtr(Vm, dict, initialDeltaT)
    );
}

} // End namespace Foam
```

- [ ] **Step 1.3: Add entries to `src/electroModels/Make/files`**

Current content of `Make/files`:
```
electroModel/electroModel.C

ecgModel/ecgModel.C
pseudoECGElectro/pseudoECGElectro.C
EikonalSolver/EikonalSolver.C
MonoDomainSolver/MonoDomainSolver.C
SingleCellSolver/SingleCellSolver.C
electroMechanicalModel/electroMechanicalModel.C


LIB = $(FOAM_MODULE_LIBBIN)/libelectroModels
```

Add after `ecgModel/ecgModel.C`:
```
purkinjeModel/purkinjeModel.C
purkinjeNetworkModel/purkinjeNetworkModel.C
```

- [ ] **Step 1.4: Verify the base class compiles (stub `purkinjeNetworkModel.C` must exist first)**

Create a temporary stub so wmake can link:

```cpp
// src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.H
#ifndef purkinjeNetworkModel_H
#define purkinjeNetworkModel_H
#include "purkinjeModel.H"
namespace Foam {
class purkinjeNetworkModel : public purkinjeModel {
public:
    TypeName("purkinjeNetworkModel");
    purkinjeNetworkModel(const volScalarField&, const dictionary&, scalar);
    void evolve(scalar, scalar, const volScalarField&, volScalarField&) override {}
};
} // End namespace Foam
#endif
```

```cpp
// src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.C
#include "purkinjeNetworkModel.H"
#include "addToRunTimeSelectionTable.H"
namespace Foam {
defineTypeNameAndDebug(purkinjeNetworkModel, 0);
addToRunTimeSelectionTable(purkinjeModel, purkinjeNetworkModel, dictionary);
purkinjeNetworkModel::purkinjeNetworkModel
(const volScalarField& Vm, const dictionary& dict, scalar dt)
: purkinjeModel(Vm, dict) {}
} // End namespace Foam
```

- [ ] **Step 1.5: Build and verify no errors**

```bash
cd /Users/simaocastro/cardiacFoamv2/src/electroModels
wmake 2>&1 | tail -20
```

Expected: `Making dependency list for source file purkinjeModel/purkinjeModel.C` and `Making dependency list for source file purkinjeNetworkModel/purkinjeNetworkModel.C` with no errors. The library `libelectroModels.so` (or `.dylib`) should be updated.

- [ ] **Step 1.6: Commit**

```bash
cd /Users/simaocastro/cardiacFoamv2
git add src/electroModels/purkinjeModel/ \
        src/electroModels/purkinjeNetworkModel/ \
        src/electroModels/Make/files
git commit -m "feat: add purkinjeModel abstract base and stub concrete class"
```

---

## Task 2: `purkinjeNetworkModel` — full header and constructor

**Files:**
- Modify: `src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.H` (replace stub)
- Modify: `src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.C` (replace stub)

This task replaces the stubs with the real class definition and a constructor that reads the graph, PVJ mapping, root stimulus, and ionic model from the dict.

- [ ] **Step 2.1: Replace `purkinjeNetworkModel.H` with full definition**

```cpp
// src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.H
#ifndef purkinjeNetworkModel_H
#define purkinjeNetworkModel_H

#include "purkinjeModel.H"
#include "ionicModel.H"

namespace Foam
{

class purkinjeNetworkModel : public purkinjeModel
{
    // ---- Graph topology (built from 'edges' list in dict) ----
    label       nNodes_;        // total node count in the Purkinje tree
    label       nEdges_;
    labelList   edgeNodeA_;     // first endpoint of each edge
    labelList   edgeNodeB_;     // second endpoint of each edge
    scalarField edgeLength_;    // L  [m]
    scalarField edgeConductance_; // sigma [S/m]

    // ---- PVJ mapping (parallel lists) ----
    labelList   pvjNodes_;      // Purkinje node indices that are PVJs
    labelList   pvjCellIDs_;    // matching 3D mesh cell indices

    // ---- Root stimulus ----
    scalar      rootStartTime_;
    scalar      rootDuration_;
    scalar      rootIntensity_; // [A/m^3] point-source applied to node 0

    // ---- Monodomain parameters for the Purkinje network ----
    scalar      chi_;           // surface-to-volume ratio [1/m]
    scalar      Cm_;            // membrane capacitance [F/m^2]

    // ---- Coupling ----
    scalar      R_pvj_;         // global PVJ resistance [Ohm.m^2]

    // ---- 1D state (size = nNodes_) ----
    scalarField Vm1D_;          // transmembrane voltage [V]
    scalarField Iion1D_;        // ionic current [V/s] (same convention as 3D Iion_)

    // ---- Ionic model (RTST, nNodes_ integration points) ----
    autoPtr<ionicModel> ionicModelPtr_;

    // ---- Private helpers ----
    void readGraph(const dictionary& dict);
    void readPVJs(const dictionary& dict);
    void readRootStimulus(const dictionary& dict);

    // No copy
    purkinjeNetworkModel(const purkinjeNetworkModel&) = delete;
    void operator=(const purkinjeNetworkModel&) = delete;

public:
    TypeName("purkinjeNetworkModel");

    purkinjeNetworkModel
    (
        const volScalarField& Vm,
        const dictionary&     dict,
        const scalar          initialDeltaT
    );

    virtual ~purkinjeNetworkModel() = default;

    virtual void evolve
    (
        scalar                t0,
        scalar                dt,
        const volScalarField& Vm,
        volScalarField&       externalStimulusCurrent
    ) override;
};

} // End namespace Foam

#endif
```

- [ ] **Step 2.2: Replace `purkinjeNetworkModel.C` with constructor implementation**

```cpp
// src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.C
#include "purkinjeNetworkModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(purkinjeNetworkModel, 0);
addToRunTimeSelectionTable(purkinjeModel, purkinjeNetworkModel, dictionary);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void purkinjeNetworkModel::readGraph(const dictionary& dict)
{
    // Each entry in 'edges' is a list of 4 scalars:
    // ( nodeA  nodeB  length  conductance )
    const List<scalarList> edgeEntries(dict.lookup("edges"));
    nEdges_ = edgeEntries.size();

    edgeNodeA_.setSize(nEdges_);
    edgeNodeB_.setSize(nEdges_);
    edgeLength_.setSize(nEdges_);
    edgeConductance_.setSize(nEdges_);

    label maxNode = 0;
    forAll(edgeEntries, eI)
    {
        const scalarList& e = edgeEntries[eI];
        if (e.size() != 4)
        {
            FatalErrorInFunction
                << "Each entry in 'edges' must have 4 values: "
                << "(nodeA nodeB length conductance). Entry " << eI
                << " has " << e.size() << " values."
                << exit(FatalError);
        }
        edgeNodeA_[eI]      = static_cast<label>(e[0]);
        edgeNodeB_[eI]      = static_cast<label>(e[1]);
        edgeLength_[eI]     = e[2];
        edgeConductance_[eI]= e[3];
        maxNode = max(maxNode, max(edgeNodeA_[eI], edgeNodeB_[eI]));
    }

    nNodes_ = maxNode + 1;

    Info<< "Purkinje network: " << nNodes_ << " nodes, "
        << nEdges_ << " edges." << nl;
}

void purkinjeNetworkModel::readPVJs(const dictionary& dict)
{
    pvjNodes_   = labelList(dict.lookup("pvjNodes"));
    pvjCellIDs_ = labelList(dict.lookup("pvjCellIDs"));

    if (pvjNodes_.size() != pvjCellIDs_.size())
    {
        FatalErrorInFunction
            << "pvjNodes and pvjCellIDs must have the same size. "
            << "pvjNodes.size()=" << pvjNodes_.size()
            << " pvjCellIDs.size()=" << pvjCellIDs_.size()
            << exit(FatalError);
    }

    Info<< "Purkinje PVJs: " << pvjNodes_.size() << " junctions." << nl;
}

void purkinjeNetworkModel::readRootStimulus(const dictionary& dict)
{
    const dictionary& rsDict = dict.subDict("rootStimulus");
    rootStartTime_ = rsDict.get<scalar>("startTime");
    rootDuration_  = rsDict.get<scalar>("duration");
    rootIntensity_ = rsDict.get<scalar>("intensity");

    Info<< "Purkinje root stimulus: start=" << rootStartTime_
        << " duration=" << rootDuration_
        << " intensity=" << rootIntensity_ << nl;
}

// * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * //

purkinjeNetworkModel::purkinjeNetworkModel
(
    const volScalarField& Vm,
    const dictionary&     dict,
    const scalar          initialDeltaT
)
:
    purkinjeModel(Vm, dict),
    nNodes_(0),
    nEdges_(0),
    edgeNodeA_(),
    edgeNodeB_(),
    edgeLength_(),
    edgeConductance_(),
    pvjNodes_(),
    pvjCellIDs_(),
    rootStartTime_(0),
    rootDuration_(0),
    rootIntensity_(0),
    chi_(dict.subDict("purkinjeNetworkModelCoeffs").get<scalar>("chi")),
    Cm_(dict.subDict("purkinjeNetworkModelCoeffs").get<scalar>("Cm")),
    R_pvj_(dict.get<scalar>("R_pvj")),
    Vm1D_(),
    Iion1D_(),
    ionicModelPtr_()
{
    readGraph(dict);
    readPVJs(dict);
    readRootStimulus(dict);

    // Initialise 1D state
    Vm1D_.setSize(nNodes_, -0.084);   // Stewart resting potential [V]
    Iion1D_.setSize(nNodes_, 0.0);

    // Construct ionic model with nNodes_ integration points
    const dictionary& coeffsDict =
        dict.subDict("purkinjeNetworkModelCoeffs");

    ionicModelPtr_ = ionicModel::New(coeffsDict, nNodes_, initialDeltaT);

    Info<< "purkinjeNetworkModel constructed with ionic model "
        << ionicModelPtr_->type() << nl << endl;
}

// * * * * * * * * * * * * * * evolve (stub) * * * * * * * * * * * * * * * * //

void purkinjeNetworkModel::evolve
(
    scalar                /*t0*/,
    scalar                /*dt*/,
    const volScalarField& /*Vm*/,
    volScalarField&       /*externalStimulusCurrent*/
)
{
    // Implemented in Task 3
}

} // End namespace Foam
```

- [ ] **Step 2.3: Build**

```bash
cd /Users/simaocastro/cardiacFoamv2/src/electroModels
wmake 2>&1 | tail -20
```

Expected: clean build. If `ionicModel::New` cannot find `Stewart`, verify that
`purkinjeNetworkModelCoeffs { ionicModel Stewart; }` is in the dict and that
`libionicModels` is linked (it already is via `Make/options`).

- [ ] **Step 2.4: Commit**

```bash
cd /Users/simaocastro/cardiacFoamv2
git add src/electroModels/purkinjeNetworkModel/
git commit -m "feat: purkinjeNetworkModel constructor reads graph, PVJs, root stimulus"
```

---

## Task 3: `purkinjeNetworkModel::evolve()` — 1D solver loop

**Files:**
- Modify: `src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.C`

Replace the stub `evolve()` with the full operator-splitting implementation.

- [ ] **Step 3.1: Replace the stub `evolve()` in `purkinjeNetworkModel.C`**

Replace:
```cpp
void purkinjeNetworkModel::evolve
(
    scalar                /*t0*/,
    scalar                /*dt*/,
    const volScalarField& /*Vm*/,
    volScalarField&       /*externalStimulusCurrent*/
)
{
    // Implemented in Task 3
}
```

With:
```cpp
void purkinjeNetworkModel::evolve
(
    scalar                t0,
    scalar                dt,
    const volScalarField& Vm3D,
    volScalarField&       externalStimulusCurrent
)
{
    // --- Step 1: Apply root stimulus at node 0 ---
    // The root stimulus is a point-source current [A/m^3] added to the
    // RHS of node 0's monodomain equation via the external current term.
    // We accumulate it into a temporary rootCurrent array and apply below.
    scalarField rootCurrent(nNodes_, 0.0);
    if (t0 >= rootStartTime_ && t0 <= (rootStartTime_ + rootDuration_))
    {
        rootCurrent[0] = rootIntensity_;
    }

    // --- Step 2: Advance ionic model for all Purkinje nodes ---
    // ionicModel::solveODE takes scalarField& directly (same convention as 3D).
    // Vm1D_ has size nNodes_; solveODE updates Iion1D_ in-place.
    ionicModelPtr_->solveODE(t0, dt, Vm1D_, Iion1D_);

    // --- Step 3: Compute 1D graph Laplacian and advance Vm1D_ explicitly ---
    // Discrete Laplacian at node i:
    //   L_i = sum_j  sigma_ij * (Vm1D_[j] - Vm1D_[i]) / L_ij^2
    // Units: [S/m * V / m^2] = [A/m^3]  (matches fvc::laplacian in 3D)
    //
    // Explicit Euler update (same form as 3D):
    //   chi*Cm * dVm/dt = L_i - chi*Cm*Iion1D_[i] + rootCurrent[i]
    //   => Vm_new[i] = Vm_old[i] + dt/(chi*Cm) * (L_i + rootCurrent[i])
    //                             - dt * Iion1D_[i]

    scalarField laplacian(nNodes_, 0.0);
    for (label eI = 0; eI < nEdges_; eI++)
    {
        const label  A     = edgeNodeA_[eI];
        const label  B     = edgeNodeB_[eI];
        const scalar L     = edgeLength_[eI];
        const scalar sigma = edgeConductance_[eI];
        const scalar flux  = sigma * (Vm1D_[B] - Vm1D_[A]) / (L * L);
        laplacian[A] += flux;
        laplacian[B] -= flux;   // symmetric: flux flows from A to B
    }

    const scalar dtOverChiCm = dt / (chi_ * Cm_);
    for (label nI = 0; nI < nNodes_; nI++)
    {
        Vm1D_[nI] += dtOverChiCm * (laplacian[nI] + rootCurrent[nI])
                   - dt * Iion1D_[nI];
    }

    // --- Step 4: Scatter coupling current to 3D externalStimulusCurrent ---
    // Vergara-Quarteroni: I_coupling = (Vm1D_pvj - Vm3D_pvj) / R_pvj
    // This is added to externalStimulusCurrent (which was already set by
    // updateExternalStimulusCurrent before this call).
    scalarField& extI = externalStimulusCurrent.primitiveFieldRef();

    forAll(pvjNodes_, k)
    {
        const label pn  = pvjNodes_[k];
        const label cId = pvjCellIDs_[k];
        const scalar Icoupling =
            (Vm1D_[pn] - Vm3D[cId]) / R_pvj_;
        extI[cId] += Icoupling;
    }

    externalStimulusCurrent.correctBoundaryConditions();
}
```

- [ ] **Step 3.2: Build**

```bash
cd /Users/simaocastro/cardiacFoamv2/src/electroModels
wmake 2>&1 | tail -20
```

Expected: clean build, no errors.

- [ ] **Step 3.3: Commit**

```bash
cd /Users/simaocastro/cardiacFoamv2
git add src/electroModels/purkinjeNetworkModel/purkinjeNetworkModel.C
git commit -m "feat: purkinjeNetworkModel::evolve — 1D graph Laplacian + coupling current"
```

---

## Task 4: Wire `purkinjeModelPtr_` into `MonoDomainSolver`

**Files:**
- Modify: `src/electroModels/MonoDomainSolver/MonoDomainSolver.H`
- Modify: `src/electroModels/MonoDomainSolver/MonoDomainSolver.C`

- [ ] **Step 4.1: Add `#include` and data member to `MonoDomainSolver.H`**

After the existing `#include "ecgModel.H"` line (around line 47):
```cpp
#include "purkinjeModel.H"
```

After the existing `autoPtr<ecgModel> ecgModelPtr_;` data member (around line 122):
```cpp
//- Optional 1D Purkinje Network model
autoPtr<purkinjeModel> purkinjeModelPtr_;
```

- [ ] **Step 4.2: Construct `purkinjeModelPtr_` in `MonoDomainSolver.C`**

In the constructor body, after the block that constructs `ecgModelPtr_`
(which ends around line 387):

```cpp
if (cardiacProperties_.found("purkinjeNetwork"))
{
    purkinjeModelPtr_ = purkinjeModel::New
    (
        Vm_,
        cardiacProperties_.subDict("purkinjeNetwork"),
        runTime.deltaTValue()
    );
}
```

- [ ] **Step 4.3: Call `purkinjeModelPtr_->evolve()` in `evolveExplicit()`**

In `evolveExplicit()`, after the call to `updateExternalStimulusCurrent` (around line 170)
and before the `ionicModelPtr_->solveODE(...)` call, add:

```cpp
// Advance 1D Purkinje network and inject coupling current
if (purkinjeModelPtr_)
{
    purkinjeModelPtr_->evolve(t0, dt, Vm_, externalStimulusCurrent_);
}
```

- [ ] **Step 4.4: Call `purkinjeModelPtr_->evolve()` in `evolveImplicit()`**

Apply the identical insertion in `evolveImplicit()`, after the call to
`updateExternalStimulusCurrent` (around line 206):

```cpp
// Advance 1D Purkinje network and inject coupling current
if (purkinjeModelPtr_)
{
    purkinjeModelPtr_->evolve(t0, dt, Vm_, externalStimulusCurrent_);
}
```

- [ ] **Step 4.5: Build**

```bash
cd /Users/simaocastro/cardiacFoamv2/src/electroModels
wmake 2>&1 | tail -20
```

Expected: clean build, `libelectroModels` updated.

- [ ] **Step 4.6: Build the solver application that links `libelectroModels`**

Locate and rebuild the main application to confirm linkage is clean:
```bash
cd /Users/simaocastro/cardiacFoamv2/applications
wmake all 2>&1 | tail -20
```

Expected: no undefined symbol errors.

- [ ] **Step 4.7: Commit**

```bash
cd /Users/simaocastro/cardiacFoamv2
git add src/electroModels/MonoDomainSolver/MonoDomainSolver.H \
        src/electroModels/MonoDomainSolver/MonoDomainSolver.C
git commit -m "feat: wire purkinjeModelPtr_ into MonoDomainSolver"
```

---

## Task 5: Integration test — minimal tutorial with Purkinje network

**Goal:** Run the existing `singleCell` or `NiedererEtAl2012` tutorial with a minimal
`purkinjeNetwork` subdict to confirm the solver advances without crashing and produces
a nonzero coupling current at the PVJ cell.

- [ ] **Step 5.1: Pick a target tutorial**

Use `tutorials/NiedererEtAl2012` — it runs a 3D monodomain solve and has a controlled
mesh. Check the mesh cell count to pick a valid `pvjCellIDs` value:

```bash
grep -r "nCells" tutorials/NiedererEtAl2012/constant/polyMesh/owner 2>/dev/null | head -3
# OR simply run:
ls tutorials/NiedererEtAl2012/constant/polyMesh/
```

The mesh has a known cell count. Pick cell index 0 as the PVJ cell for the test
(always valid regardless of mesh size).

- [ ] **Step 5.2: Add `purkinjeNetwork` subdict to the tutorial's `electroProperties`**

Edit `tutorials/NiedererEtAl2012/constant/electroProperties`.
Add this block at the end, before the closing `}`:

```cpp
purkinjeNetwork
{
    purkinjeModel    purkinjeNetworkModel;

    R_pvj            500.0;   // [Ohm.m^2]

    // Minimal 2-node, 1-edge tree: root=0, terminal=1
    edges
    (
        ( 0  1  0.001  0.003 )
    );

    pvjNodes   ( 1 );
    pvjCellIDs ( 0 );    // cell 0 is always valid

    rootStimulus
    {
        startTime    0.0;
        duration     0.002;
        intensity    80000.0;
    }

    purkinjeNetworkModelCoeffs
    {
        ionicModel    Stewart;

        chi            1400.0;    // [1/m]  typical Purkinje value
        Cm             1.0e-2;    // [F/m^2]
    }
}
```

- [ ] **Step 5.3: Run the tutorial and verify no crash**

```bash
cd tutorials/NiedererEtAl2012
./Allclean
./Allrun 2>&1 | head -60
```

Expected output contains:
```
Selecting purkinjeModel purkinjeNetworkModel
Purkinje network: 2 nodes, 1 edges.
Purkinje PVJs: 1 junctions.
Purkinje root stimulus: start=0 duration=0.002 intensity=80000
purkinjeNetworkModel constructed with ionic model Stewart
```

And the solver runs to completion without segfault or `FatalError`.

- [ ] **Step 5.4: Verify nonzero coupling current is produced**

After the run, check that `externalStimulusCurrent` at cell 0 is nonzero at `t=0.001`
(within the root stimulus window). The fastest way:

```bash
# Check the field at the first output time
grep -A 5 "internalField" tutorials/NiedererEtAl2012/0.001/externalStimulusCurrent 2>/dev/null
```

If `externalStimulusCurrent` is not written, add `AUTO_WRITE` flag temporarily or print
`externalStimulusCurrent[0]` with an `Info<<` statement in `evolveExplicit()` for
one timestep. The value should differ from the case without Purkinje.

- [ ] **Step 5.5: Restore tutorial to clean state (keep `electroProperties` with Purkinje)**

The modified `electroProperties` (with `purkinjeNetwork`) is the new reference config.
Copy it into `referenceTest/` so regression tests lock it:

```bash
cp tutorials/NiedererEtAl2012/constant/electroProperties \
   tutorials/NiedererEtAl2012/referenceTest/constant/electroProperties
```

- [ ] **Step 5.6: Commit**

```bash
cd /Users/simaocastro/cardiacFoamv2
git add tutorials/NiedererEtAl2012/constant/electroProperties \
        tutorials/NiedererEtAl2012/referenceTest/constant/electroProperties
git commit -m "test: add minimal purkinjeNetwork to NiedererEtAl2012 tutorial"
```

---

## Self-Review Notes

- **Spec coverage:** All 5 design sections are covered: class structure (Task 1), dict format (Task 2), internal state (Task 2), call sequence (Tasks 3–4), files (all tasks).
- **`ionicModel::New` dict key:** `ionicModel::New` reads `ionicModel` from the passed dict. The `purkinjeNetworkModelCoeffs` subdict contains `ionicModel Stewart;` — this matches.
- **`Vm3D[cId]` access:** `volScalarField` supports `operator[]` on `label` returning the cell value — this is valid OpenFOAM.
- **`externalStimulusCurrent.primitiveFieldRef()`:** Returns non-const `scalarField&` to internal data — correct for in-place modification before `correctBoundaryConditions()`.
- **Laplacian sign:** The loop adds `+flux` to A and `-flux` to B. For `flux = sigma*(Vm_B - Vm_A)/L^2`, this means node A gains positive current when `Vm_B > Vm_A` — consistent with diffusion from high to low potential.
- **Stewart resting potential:** `Vm1D_` initialised at `-0.084` V to match the `Vm_` initial condition in `MonoDomainSolver`.
