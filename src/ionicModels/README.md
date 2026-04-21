# ionicModels library architecture

This directory provides runtime-selectable ionic cell models (`libionicModels`) used by
myocardium workflows, Purkinje/conduction workflows, and `singleCellSolver`.

Each model wraps generated ODE code in a shared `Foam::ionicModel` interface.

## Directory structure

```text
src/ionicModels/
├── ionicModel/            # Base class, selector utilities, runtime selection
├── monodomainFDAManufactured/  # Manufactured-solution verification model
├── bidomainFDAManufactured/    # Manufactured-solution verification model
├── AlievPanfilov/
├── BuenoOrovio/
├── Courtemanche/
├── Fabbri/
├── Gaur/
├── Grandi/
├── ORd/
├── Stewart/
├── TNNP/
├── ToRORd_dynCl/
├── Trovato/
├── Make/
└── README.md
```

## Core class: `Foam::ionicModel`

Defined in `ionicModel/ionicModel.H` and implemented in `ionicModel/ionicModel.C`.

### Main responsibilities

- Inherits `ODESystem` (OpenFOAM ODE API).
- Implements runtime selection via `ionicModel::New(...)`.
- Stores per-integration-point ODE step sizes (`step_`).
- Stores tissue/dimension selector flag (`tissue_`).
- Supports single-cell mode (`solveVmWithinODESolver_`).
- Stores stimulation protocol (`StimulusProtocol`) loaded from dictionary.
- Implements optional generic export/write/debug APIs through `ionicModelIO`.
- Implements `ElectromechanicalSignalProvider` interface for electromechanics signal exchange.

### Key virtual interface

Derived models implement:

- `solveODE(...)`
- `derivatives(...)`
- `nEqns() const`

Optional overrides:

- `supportedTissueTypes()`
- `supportedDimensions()`
- `sweepCurrent(...)`
- coupling-signal methods (`hasSignal`, `signal`)

### Coupling signals in the base class

`ionicModel` now provides default coupling signal exposure from model metadata:

- `Vm`: available when a model provides `ioVmTransform()` or has a voltage state
  (e.g. `membrane_V`, `V`, `Vm`).
- `Cai`: available when a state matching intracellular calcium naming exists
  (e.g. `Ca_i`, `Cai`, `calcium_Cai`).

This means most detailed ionic models expose `Vm`/`Cai` without per-model
signal boilerplate.

## Tissue vs dimension selection

`ionicModel/ionicSelector` centralizes dictionary interpretation:

- Normal models use `tissue` entry.
- Manufactured/verification models can use `dimension` entry via model-side selection.

This keeps selection logic consistent across all ionic models.

## I/O and export architecture

The base class can write/export without model-specific code when metadata hooks are provided
(`ioStateNames`, `ioAlgebraicNames`, `ioStatesPtr`, etc.).

Common behaviors:

- Filter exported/debug variable lists.
- Write full or selected headers.
- Import solver-owned volumetric fields back into ionic state storage when needed.
- Export selected variables into `volScalarField` lists.
- Support relaxed variable name compatibility for Vm/rates through `ionicVariableCompatibility`.

## Compiled ionic models

Current `Make/files` entries:

- `AlievPanfilov`
- `BuenoOrovio`
- `Courtemanche`
- `Fabbri`
- `Gaur`
- `Grandi`
- `ORd`
- `Stewart`
- `TNNP`
- `ToRORd_dynCl`
- `Trovato`
- `monodomainFDAManufactured`
- `bidomainFDAManufactured`

## Build target

`Make/files` builds into:

- `$(FOAM_USER_LIBBIN)/libionicModels`

## Adding a new ionic model

1. Add model folder with generated equations and wrapper `.H/.C`.
2. Derive from `Foam::ionicModel` and implement required methods.
3. Register runtime type in the model `.C` file.
4. Add the `.C` file to `src/ionicModels/Make/files`.
5. Provide metadata hooks if generic export/write behavior is desired.
