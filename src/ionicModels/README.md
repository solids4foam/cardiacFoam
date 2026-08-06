# ionicModels library architecture

This directory provides runtime-selectable ionic cell models (`libionicModels`) used by
myocardium workflows, Purkinje/conduction workflows, and `singleCellSolver`.

Each model wraps generated ODE code in a shared `Foam::ionicModel` interface.

## Directory structure

```text
src/ionicModels/
├── ionicModel/                              # Core classes and utilities
│   ├── ionicModel.H/.C                      # Base class; virtual interface
│   ├── configuredIonicModel.H               # Scalar model heterogeneity layer
│   ├── configuredBatchedIonicModel.H        # Batched model heterogeneity layer
│   ├── ionicHeterogeneity.H/.C              # Transmural/regional weight calculations
│   ├── ionicHeterogeneityOrchestrator.H/.C # Heterogeneity mode dispatch
│   ├── ionicSelector.H/.C                   # Tissue/dimension dictionary parsing
│   └── batchedIonicModel.H/.C               # GPU-oriented base
├── verificationModels/
│   ├── monodomainFDAManufactured/   # Manufactured-solution verification model
│   ├── bidomainFDAManufactured/     # Manufactured-solution verification model
│   └── bathBidomainFDAManufactured/ # Manufactured-solution verification model (bath bidomain)
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

Most ionic models expose `Vm`/`Cai` automatically via base-class metadata.

## Tissue vs dimension selection

`ionicModel/ionicSelector` centralizes dictionary interpretation:

- Normal models use `tissue` entry.
- Manufactured/verification models can use `dimension` entry via model-side selection.

This keeps selection logic consistent across all ionic models.

## I/O and export architecture

The base class can write/export without model-specific code when metadata hooks are provided
(`ioStateNames`, `ioAlgebraicNames`, `ioStatesPtr`, etc.).

Common behaviors:

- Filtering of exported/debug variable lists
- Full or selected header writing
- Import of solver-owned volumetric fields back into ionic state storage when needed
- Export of selected variables into `volScalarField` lists
- Relaxed variable name compatibility for Vm/rates through `ionicVariableCompatibility`

## Runtime model contract

The values below are the exact `ionicModel` dictionary selectors registered by
the sources in `Make/files`. All are built in both full and lightweight modes.
CUDA kernels are added only when `CARDIAC_ENABLE_CUDA` is set; compact-batched
models otherwise use their maintained host backend.

| Purpose | Runtime name | Backend | Source boundary | Known equivalence limitation |
|---|---|---|---|---|
| Cell model | `AlievPanfilov` | scalar CPU | maintained wrapper; generated ODE/Names headers | Not numerically identical to compact-batched integration by construction |
| Cell model | `BuenoOrovio` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `Courtemanche` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `Fabbri` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `Gaur` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `Grandi` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `PerisYague` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `Stewart` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `TNNP` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `ToRORd_dynCl` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `Trovato` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `TWorld` | scalar CPU | maintained wrapper; generated ODE/Names headers | As above |
| Cell model | `AlievPanfilovcompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `AlievPanfilovBatched` is not a runtime selector |
| Cell model | `BuenoOroviocompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `BuenoOrovioBatched` is not a runtime selector |
| Cell model | `CourtemanchecompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `CourtemancheBatched` is not a runtime selector |
| Cell model | `FabbricompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `FabbriBatched` is not a runtime selector |
| Cell model | `GaurcompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `GaurBatched` is not a runtime selector |
| Cell model | `GrandicompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `GrandiBatched` is not a runtime selector |
| Cell model | `PerisYaguecompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `PerisYagueBatched` is not a runtime selector |
| Cell model | `StewartcompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `StewartBatched` is not a runtime selector |
| Cell model | `TNNPcompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `TNNPBatched` is not a runtime selector |
| Cell model | `ToRORd_dynClcompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `ToRORd_dynClBatched` is not a runtime selector |
| Cell model | `TrovatocompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `TrovatoBatched` is not a runtime selector |
| Cell model | `TWorldcompactBatched` | SoA host; optional CUDA | maintained wrapper/backend; generated batch equations | `TWorldBatched` is not a runtime selector |
| Monodomain MMS forcing | `monodomainFDAManufactured` | scalar CPU | maintained verification implementation | Verification-only; not a physiological model |
| Bidomain MMS forcing | `bidomainFDAManufactured` | scalar CPU | maintained verification implementation | Verification-only; not a physiological model |
| Bath-bidomain MMS forcing | `bathBidomainFDAManufactured` | scalar CPU | maintained verification implementation | Verification-only; not a physiological model |

### Generated and maintained boundaries

- Wrapper `.C`/`.H` files, runtime registration, dictionary handling, host/CUDA
  dispatch, and heterogeneity orchestration are maintained project code.
- Year-labelled equation, Names, and Batch headers are generated model code.
  Change their generator/mapping contract and regenerate rather than hand-editing
  equations or metadata in isolation.
- Files listed by `Make/files-gpu` are maintained CUDA integration kernels. Their
  presence does not create additional runtime names.
- Scalar and compact-batched variants implement the same named model family but
  use different layouts and integration paths. Exact trajectory equivalence is
  not guaranteed; compare with explicit tolerances before changing references.

## Build target

`Make/files` builds into:

- `$(FOAM_USER_LIBBIN)/libionicModels`

## Tissue heterogeneity

Spatial heterogeneity of ionic properties is supported through the optional
`ionicHeterogeneity` dictionary block. All non-batched scalar models inherit
from `configuredIonicModel` and support region-based heterogeneity. Batched
models support heterogeneity only if their `supportedTissueTypes()` includes
all three anatomical tissue types (endocardialCells, mCells, epicardialCells).

### Supported models

**Scalar CPU models with native tissue baselines in the generated model:**

- `BuenoOrovio`
- `TNNP`
- `ToRORd_dynCl`
- `TWorld`

These models branch on tissue selection in their built-in constants and/or
initial states before any override is applied.

**Scalar CPU models with override-driven region support from one baseline:**

- `AlievPanfilov`
- `Courtemanche`
- `Fabbri`
- `Gaur`
- `Grandi`
- `PerisYague`
- `Stewart`
- `Trovato`

These models accept region labels and spatial heterogeneity, but the generated
core does not provide distinct endo/M/epi baselines. Region differences come
from `ionicConstantOverrides` and optional apex-to-base scaling applied on top
of the model baseline.

**Batched/GPU models (heterogeneity support where tissue types permit):**

- `BuenoOroviocompactBatched` (supports all three tissue types)
- `TNNPcompactBatched` (supports all three tissue types)
- `ToRORd_dynClcompactBatched` (supports all three tissue types)
- `TWorldcompactBatched` (supports all three tissue types)
- Other batched models (`AlievPanfilovcompactBatched`, `CourtemanchecompactBatched`, `FabbricompactBatched`, `GaurcompactBatched`, `GrandicompactBatched`, `PerisYaguecompactBatched`, `StewartcompactBatched`, `TrovatocompactBatched`) support only myocyte tissue and do not support heterogeneity

### Configuration

The `ionicHeterogeneity` block is nested within the model coefficients
(e.g., `monodomainSolverCoeffs`). The mode selects how regions are defined
and applied.

#### Mode: `transmuralBands` (default)

Classic transmural heterogeneity: partitions the wall by normalized transmural distance.

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

| Key | Default | Accepted values |
|-----|---------|-----------------|
| `field` | `t` | Any scalar field name (0=endo, 1=epi) |
| `endoMInterface` | `0.3` | 0 < value < `mEpiInterface` |
| `mEpiInterface` | `0.7` | `endoMInterface` < value < 1 |
| `transitionWidth` | `0.1` | ≥ 0; must not cause overlapping bands in blend mode |
| `transitionMode` | `blend` | `blend`, `hard` |
| `smoothing` | `smoothstep` | `smoothstep` |

**Validation:** 0 < `endoMInterface` < `mEpiInterface` < 1, and in blend mode both interfaces must respect transition width constraints.

#### Mode: `namedRegions`

Open dictionary of arbitrary named regions, sorted by field range. Each region carries a baseline tissue and optional overrides scope name.

```
ionicHeterogeneity
{
    field             t;                  // Name of field (0 at one extreme, 1 at other)
    mode              namedRegions;
    transitionWidth   0.1;                // Smooth transition band width
    transitionMode    blend;              // Transition type: blend or hard
    smoothing         smoothstep;         // Smoothing function: smoothstep

    regions
    {
        subendo
        {
            range       (0.0 0.3);
            baseline    endocardialCells;  // Default tissue type for constants
        }
        midwall
        {
            range       (0.3 0.7);
            baseline    mCells;
        }
        subepi
        {
            range       (0.7 1.0);
            baseline    epicardialCells;
        }
    }
}
```

Each region must declare `range (min max)`. Optional `baseline` defaults to the region name if it is one of `{endocardialCells, mCells, epicardialCells}`, otherwise defaults to `myocyte`. Regions must be contiguous (no gaps/overlaps) and tile [0, 1] exactly.

#### Mode: `cellZoneRegions`

Regions defined by mesh cell zone membership. Each cell is assigned exactly one region (no blending).

```
ionicHeterogeneity
{
    mode cellZoneRegions;

    regions
    {
        infarcted
        {
            cellZone    scar_zone;
            baseline    myocyte;            // Optional; defaults per baseline rule
        }
        border_zone
        {
            cellZone    border_zone_cells;
            baseline    epicardialCells;
        }
    }
}
```

Each region must declare `cellZone <meshCellZoneName>`. Optional `baseline` defaults per the standard rule. No transition width or smoothing; assignment is hard (binary per cell).

#### Apex-to-base exponential scaling

Optional `apexBaseBands` sub-dictionary (applies with any mode) scales selected constants along an apex-to-base gradient:

```
ionicHeterogeneity
{
    mode namedRegions;  // or transmuralBands, or cellZoneRegions
    regions { ... }

    apexBaseBands
    {
        beta        3.0;        // Exponential power (default)
        scalingMin  0.2;        // Scale at apex (default)
        scalingMax  5.0;        // Scale at base (default)
        variables   (G_K1 G_Na);  // Constant names to scale
        field       d;          // Apex-to-base distance field
    }
}
```

### Tissue baseline and override scoping

Each region (in any mode) specifies an optional `baseline` tissue keyword:

- Explicit: one of `epicardialCells`, `mCells`, `endocardialCells`, `myocyte`
- Default rule: if the region name itself is one of those three anatomical names, use it; otherwise `myocyte`

The baseline determines which `constantsForTissue()` and
`initialStatesForTissue()` values are used as the region's foundation.
Overrides are then applied by scope name (see below). For models without
native endo/M/epi branches, the baseline may still be `epicardialCells`,
`mCells`, or `endocardialCells`, but those labels only select the override
path; they do not imply distinct built-in CellML tissue families.

### Constant override scoping (`ionicConstantOverrides`)

Constant overrides are stored in a top-level `ionicConstantOverrides` block with named sub-dictionaries representing scopes:

```
monodomainSolverCoeffs
{
    ionicModel      BuenoOrovio;

    ionicHeterogeneity { ... }

    ionicConstantOverrides
    {
        global           // Global scope (applied first)
        {
            G_K1   0.123;
            G_Na   0.456;
        }
        epicardialCells  // Tissue-name scope
        {
            G_K1   0.130;
        }
        myocyte          // Fallback tissue name
        {
            G_Na   0.450;
        }
        scar_region      // Named region scope (if ionicHeterogeneity uses named regions)
        {
            G_K1   0.050;
        }
    }
}
```

Scopes are applied in order:

1. `global` (always applied first)
2. Tissue scope matching the region's `baseline` (endocardialCells/mCells/epicardialCells/myocyte)
3. Named region scope (for namedRegions/cellZoneRegions mode only)

Scopes 2 and 3 are optional. If a constant appears in multiple scopes, later scopes override earlier ones.

## Adding a new ionic model

1. Add model folder with generated equations and wrapper `.H/.C`.

2. **Choose the base class:**
   - For a production model that should support heterogeneity: derive from `Foam::configuredIonicModel` (for scalar CPU) or `Foam::configuredBatchedIonicModel` (for batched GPU). These inherit from `ionicModel` and provide `HETEROGENEOUS_CONSTANTS_` storage and forwarding overrides for free.
   - For a verification/manufactured model or one that explicitly does not support heterogeneity: derive directly from `Foam::ionicModel`. If heterogeneity is requested at runtime, it will fatal with a clear error.

3. Implement required virtual methods:
   - `solveODE(...)` / `evaluateState(...)`
   - `derivatives(...)`
   - `nEqns() const`
   - `constantsForTissue(tissueFlag)`
   - `initialStatesForTissue(tissueFlag)`

4. Override metadata hooks for export (optional):
   - `ioVmTransform()`, `ioStateNames()`, `ioStatesPtr()`, `ioConstantNames()`, etc.
   - See existing models (e.g., BuenoOrovio.H) for the pattern.

5. Optionally override `supportedTissueTypes()` and `supportedDimensions()` to advertise capabilities.
   - For `configuredIonicModel` / `configuredBatchedIonicModel` to enable heterogeneity at runtime, `supportedTissueTypes()` must include the anatomical tissue names the model supports.
   - Batched models are only enabled for heterogeneity if they return all three (endocardialCells, mCells, epicardialCells).

6. Per-cell initial-state blending is handled differently depending on the base class:
   - `configuredBatchedIonicModel` already applies `HETEROGENEOUS_INITIAL_STATES_` generically in its shared `configureIonicHeterogeneity()` — no batched model needs to override this itself (`ToRORd_dynClBatched` does not).
   - `configuredIonicModel` does *not* do this generically. If a scalar model needs per-cell initial-state blending, it must override `configureIonicHeterogeneity()` itself and apply `HETEROGENEOUS_INITIAL_STATES_` manually, as `ToRORd_dynCl` does.

7. Register runtime type via `OverrideTypeName("ModelName")` in the class definition.

8. Add the `.C` file to `src/ionicModels/Make/files`.

9. Ensure `constantsForTissue()` and `initialStatesForTissue()` fill the returned fields correctly for each tissue flag (from `ionicModel`'s tissue selector).
