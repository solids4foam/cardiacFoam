# ionicModels library architecture

This directory provides runtime-selectable ionic cell models (`libionicModels`) used by
myocardium workflows, Purkinje/conduction workflows, and `singleCellSolver`.

Each model wraps generated ODE code in a shared `Foam::ionicModel` interface.

## Directory structure

```text
src/ionicModels/
├── ionicModel/                  # Base class, selector, batched/GPU support headers
├── monodomainFDAManufactured/   # Manufactured-solution verification model
├── bidomainFDAManufactured/     # Manufactured-solution verification model
├── bathBidomainFDAManufactured/ # Manufactured-solution verification model (bath bidomain)
├── AlievPanfilov/
├── AlievPanfilovBatched/
├── BuenoOrovio/
├── BuenoOrovioBatched/
├── Courtemanche/
├── CourtemancheBatched/
├── Fabbri/
├── FabbriBatched/
├── Gaur/
├── GaurBatched/
├── Grandi/
├── GrandiBatched/
├── PerisYague/
├── PerisYagueBatched/
├── Stewart/
├── StewartBatched/
├── TNNP/
├── TNNPBatched/
├── ToRORd_dynCl/
├── ToRORd_dynClBatched/
├── Trovato/
├── TrovatoBatched/
├── TWorld/
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

**Scalar models:**

- `AlievPanfilov`
- `BuenoOrovio`
- `Courtemanche`
- `Fabbri`
- `Gaur`
- `Grandi`
- `PerisYague`
- `Stewart`
- `TNNP`
- `ToRORd_dynCl`
- `Trovato`
- `TWorld`

**Batched (SoA) models:**

- `AlievPanfilovBatched`
- `BuenoOrovioBatched`
- `CourtemancheBatched`
- `FabbriBatched`
- `GaurBatched`
- `GrandiBatched`
- `PerisYagueBatched`
- `StewartBatched`
- `TNNPBatched`
- `ToRORd_dynClBatched`
- `TrovatoBatched`
- `TWorldBatched`

**Verification models:**

- `monodomainFDAManufactured`
- `bidomainFDAManufactured`
- `bathBidomainFDAManufactured`

## Build target

`Make/files` builds into:

- `$(FOAM_USER_LIBBIN)/libionicModels`

## Tissue heterogeneity

Transmural heterogeneity of ionic properties is supported through the optional
`ionicHeterogeneity` dictionary block, which allows spatial variation of cellular
phenotypes (endo, mid-myocardial, epi) across the wall thickness.

### Supported models

**Scalar CPU models:**

- `BuenoOrovio` (only)

**Batched/GPU models:**

- `BuenoOrovioBatched`
- `TNNPBatched`
- `TWorldBatched`
- `ToRORd_dynClBatched`

Note: Other scalar models (TNNP, TWorld, ToRORd_dynCl) support tissue-dependent
constant overrides at initialization but do not implement spatial heterogeneity.

### Configuration

The `ionicHeterogeneity` block is nested within the model coefficients
(e.g., `monodomainSolverCoeffs`):

```
ionicHeterogeneity
{
    field             t;                  // Name of transmural distance field
    mode              transmuralBands;    // Heterogeneity mode
    endoMInterface    0.3;                // Endo-to-M-cell interface
    mEpiInterface     0.7;                // M-cell-to-epi interface
    transitionWidth   0.1;                // Smooth transition band width
    transitionMode    blend;              // Transition type: blend or hard
    smoothing         smoothstep;         // Smoothing function: smoothstep
}
```

**Dictionary keys:**

| Key | Meaning | Default | Accepted values |
|-----|---------|---------|-----------------|
| `field` | Name of the transmural distance field (0 at endo, 1 at epi) | `t` | Any field name |
| `mode` | Heterogeneity application mode | `transmuralBands` | `transmuralBands` |
| `endoMInterface` | Transmural position of endo/M-cell boundary (normalized [0, 1]) | `0.3` | 0 < value < `mEpiInterface` |
| `mEpiInterface` | Transmural position of M-cell/epi boundary (normalized [0, 1]) | `0.7` | `endoMInterface` < value < 1 |
| `transitionWidth` | Width of smooth transition region | `0.1` | ≥ 0; must not cause overlapping bands in blend mode |
| `transitionMode` | Hard or smooth transitions between bands | `blend` | `blend`, `hard` |
| `smoothing` | Smoothing function applied to transition zones | `smoothstep` | `smoothstep` |

**Validation rules:**

- 0 < `endoMInterface` < `mEpiInterface` < 1
- `transitionWidth` ≥ 0
- In `blend` mode: `endoMInterface + transitionWidth` ≤ `mEpiInterface` and `mEpiInterface + transitionWidth` ≤ 1

### Tissue types

The base selection uses the `tissue` entry (outside `ionicHeterogeneity`):

- `endocardialCells`
- `mCells`
- `epicardialCells`
- `myocyte` (default if not specified)

BuenoOrovio supports all three tissue types; batched models adapt heterogeneity
weights per tissue class.

## Adding a new ionic model

1. Add model folder with generated equations and wrapper `.H/.C`.
2. Derive from `Foam::ionicModel` and implement required methods.
3. Register runtime type in the model `.C` file.
4. Add the `.C` file to `src/ionicModels/Make/files`.
5. Provide metadata hooks if generic export/write behavior is desired.
