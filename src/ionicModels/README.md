# ionicModels library

The `ionicModels` library provides a run-time selectable family of cardiac
 electrophysiology (cell / ionic) models for use with **electroActivationFoam**.
 It wraps CellML-generated ODE systems in a common C++ interface
 (`Foam::ionicModel`) and connects them to monodomain/bidomain
 reaction–diffusion solvers.

The goals of the library are:

- To keep each ionic model modular and close to its source publication,
- To expose a consistent API for:
  - Evaluating ionic current densities,
  - Advancing state variables in time,
  - Exporting states, algebraic variables, and rates,
  - Performing manufactured-solution verification where available.

---

## Directory structure

```
ionicModels/
├── BuenoOrovio/        # Bueno-Orovio et al. 2008 model
├── Courtemanche/       # Courtemanche et al. 1998 model
├── Gaur/               # Gaur et al. 2021 model
├── TNNP/               # Ten Tusscher–Noble–Noble–Panfilov 2004 model
├── tmanufacturedFDA/   # Manufactured-solution / FDA verification model
├── ionicModel/         # Base class, selectors, utility functions
├── genericWriter/      # Common I/O for single-cell and 3D output
├── Make/               # Build system files
└── README.md
```

---

## Core architecture

### `Foam::ionicModel` (base class)

Located in `ionicModel/ionicModel.H`, this abstract base class defines the API
that all ionic models must implement. It:

- Derives from `Foam::ODESystem`,
- Owns an OpenFOAM `ODESolver`,
- Stores:
  - The model dictionary (`dict_`),
  - A per-integration-point time step array (`step_`),
  - A tissue-type flag (`tissue_`),
  - An optional flag enabling the solution of `Vm` inside the ODE model
    (for single-cell simulations).

The class provides:

- **Pure virtual functions** that derived models must implement:
  - `solveODE(...)`
  - `calculateCurrent(...)`
  - `derivatives(...)`
  - `nEqns() const`
- **Run-time selection:**
  - `TypeName("ionicModel")`
  - `declareRunTimeSelectionTable(...)`
  - `autoPtr<ionicModel> ionicModel::New(...)`
- **Utility functions** for copying and snapshotting internal model states:
  - `copyInternalToExternal(...)`
  - `saveStateSnapshot(...)`
  - `restoreStateSnapshot(...)`
- **Optional features:**
  - `supportedTissueTypes()`
  - `supportedDimensions()`
  - `hasManufacturedSolution()`
  - `writeHeader(...)`, `write(...)`
  - `exportStates(...)`
  - `debugPrintFields(...)`
  - `exportedFieldNames()`, `debugPrintedNames()`

All concrete ionic models inherit from this class and override the appropriate
behaviour.

---

## Available ionic models

### BuenoOrovio (2008)

```
BuenoOrovio/
├── BuenoOrovio_2008.H        # Generated ODE system
├── BuenoOrovio_2008Names.H   # Enums and variable names
├── BuenoOrovio.H, BuenoOrovio.C
├── fieldInit.H
└── heaviside.H
```

### Courtemanche (1998)

```
Courtemanche/
├── Courtemanche_1998.H
├── Courtemanche_1998Names.H
├── Courtemanche.H
└── Courtemanche.C
```

### Gaur (2021)

```
Gaur/
├── Gaur_2021.H
├── Gaur_2021Names.H
├── Gaur.H
└── Gaur.C
```

### TNNP (2004)

```
TNNP/
├── tentusscher_noble_noble_panfilov_2004.H
├── TNNP_2004.H
├── TNNP_2004Names.H
├── TNNP.C
└── TNNP.H
```

### Manufactured FDA model

```
tmanufacturedFDA/
├── tmanufacturedFDA_2014.H
├── tmanufacturedFDA_2014Names.H
├── tmanufacturedFDA.H
├── tmanufacturedFDA.C
└── tmanufacturedFields.H
```

---

## I/O utilities: `ionicModelIO`

Located in `genericWriter/`.

This module centralises:

- header writing for single-cell CSV/TSV outputs,
- writing state/algebraic/rate traces,
- exporting internal model data as `volScalarField` objects,
- generating lists of variable names for exporting and debugging.

---

## Usage inside a PDE solver

1. Allocate:
   - `volScalarField Vm;`
   - `scalarField Iion;`
   - `Field<Field<scalar>> states;`

2. On each timestep:
   - `model.calculateCurrent(t0, dt, Vm, Iion, states);`
   - or `model.solveODE(t0, dt, Vm, Iion, states);`

3. Optionally export state fields.

---

## Adding a new ionic model

Steps:
1. Generate ODE code (usually from CellML).
2. Create new directory with wrapper `.H/.C`.
3. Implement overrides of `solveODE`, `calculateCurrent`, `derivatives`, `nEqns`.
4. Register in run‑time selection table.

---

## Manufactured solution support

The `tmanufacturedFDA` model demonstrates how to:

- override `hasManufacturedSolution()`,
- initialise analytic fields,
- export error fields for convergence testing.

---

## Building

```bash
wmake libso .
```
