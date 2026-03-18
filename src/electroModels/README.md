# electroModels library architecture

This directory provides runtime-selectable electrophysiology models used by `cardiacFoam`.
The library is built as `libelectroModels`.

## Directory structure

```text
src/electroModels/
├── electroModel/               # Base class and runtime selection
├── monoDomainElectro/          # Full monodomain PDE-ODE model
├── singleCellElectro/          # Single-cell ODE-only driver
├── eikonalDiffusionElectro/    # Reduced-order activation-time model
├── electroModelECG/            # Electro + ECG wrapper model
├── electroMechanicalModel/     # Electro-mechanics wrapper model
├── Make/
└── README.md
```

## Core design

### 1) `Foam::electroModel` (base class)

Defined in `electroModel/electroModel.H`.

Responsibilities:

- Derives from `physicsModel` and `IOdictionary`.
- Reads `constant/electroProperties`.
- Exposes common mesh/time access and the `evolve()` interface.
- Owns runtime selection table (`electroModel::New(...)`).
- Parses `solutionAlgorithm` enum (`implicit` / `explicit`).

Selection key in `electroProperties`:

```cpp
electroModel monoDomainElectro;
```

### 2) `monoDomainElectro`

Tissue-scale PDE-ODE model:

- Solves transmembrane voltage `Vm` on the mesh.
- Uses runtime-selected ionic model per cell (`ionicModel::New(...)`).
- Creates an optional ionic-model pre/post processor keyed by ionic-model type
  for mesh-dependent setup/verification.
- The shared pre/post processor abstractions live in
  `src/modelPrePostProcessors/`.
- Maintains external state storage (`Field<Field<scalar>> states_`).
- Supports both explicit and implicit stepping paths.
- Tracks `activationTime` when `Vm` crosses threshold.
- Exports selected ionic variables to volumetric fields.

Stimulus architecture:

- PDE stimulus is loaded from `monodomainStimulus` dictionary entries.
- Ionic-model S1/S2 stimulus can exist for single-cell workflows.
- Guard logic prevents accidental double stimulation in monodomain runs unless explicitly allowed (`allowIonicStimulusInMonodomain`).

### 3) `singleCellElectro`

Single integration-point workflow:

- Instantiates ionic model with one integration point.
- Enables `solveVmWithinODESolver=true` so voltage is advanced inside ODEs.
- Writes traces to `postProcessing/<ionicModel>_<tissue>_<stimulus>.txt`.
- Uses shared `ionicModelIO` helpers for headers and row writes.

### 4) `eikonalDiffusionElectro`

Reduced-order activation model:

- Solves steady eikonal-diffusion formulation for activation time `psi`.
- Uses anisotropic conductivity tensor and optional stabilized form (`eikonalAdvectionDiffusionApproach`).
- Applies stimulus through geometric box selection in the mesh.

### 5) `electroModelECG`

Wrapper workflow for cases that solve electrophysiology and compute ECG outputs
in the same run:

- Owns an inner electro model plus a runtime-selected `ecgModel`.
- Reads the `ECG` sub-dictionary from `constant/electroProperties`.
- Delegates time stepping to the electro model and ECG post-processing to the
  ECG stack in `src/ecgModels/`.

### 6) `electroMechanicalModel`

Wrapper workflow for coupled electro-mechanics runs:

- Owns an electro model plus a runtime-selected active-tension model.
- Passes coupling signals from the ionic model stack into the active-tension
  stack.
- Uses the same runtime-selection pattern as the other electro models.

## Runtime chain with `cardiacFoam`

1. `cardiacFoam` creates `physicsModel`.
2. Electro `physicsModel` creates `electroModel` from `electroProperties`.
3. `monoDomainElectro` and `singleCellElectro` create `ionicModel` from the
   selected electro-model coeff dictionary in `electroProperties`.

## Build target

`Make/files` builds:

- `electroModel/electroModel.C`
- `eikonalDiffusionElectro/eikonalDiffusionElectro.C`
- `monoDomainElectro/monoDomainElectro.C`
- `singleCellElectro/singleCellElectro.C`
- `electroModelECG/electroModelECG.C`
- `electroMechanicalModel/electroMechanicalModel.C`

into:

- `$(FOAM_MODULE_LIBBIN)/libelectroModels`
