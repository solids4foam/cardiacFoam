# electroModels library architecture

This directory provides runtime-selectable electrophysiology models used by `cardiacFoam`.
The library is built as `libelectroModels`.

## Directory structure

```text
src/electroModels/
├── core/
│   └── electroModel/           # Core orchestration and base interfaces
├── domains/
│   ├── conductionSystemDomain/ # Pre-primary auxiliary domain wrapper
│   ├── ecgDomain/              # Post-primary ECG contracts + solver interface
│   └── reactionDiffusionMyocardiumDomain/ # Shared monodomain/bidomain tissue domain
├── ecgModels/                  # Concrete ECG domain model implementations
├── conductionSystemModels/     # Purkinje / graph-backed conduction models
├── solvers/
│   ├── monoDomainSolver/      # Full monodomain PDE-ODE model
│   ├── singleCellSolver/      # Single-cell ODE-only driver
│   └── eikonalSolver/         # Reduced-order activation-time model
├── wrappers/
│   └── electroMechanicalModel/ # Electro-mechanics wrapper model
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
electroModel MonoDomainSolver;
```

### 2) `MonoDomainSolver`

Tissue-scale PDE-ODE model:

- Solves transmembrane voltage `Vm` on the mesh.
- Uses runtime-selected ionic model per cell (`ionicModel::New(...)`).
- Creates an optional ionic-model pre/post processor keyed by ionic-model type
  for mesh-dependent setup/verification.
- The shared pre/post processor abstractions live in
  `src/modelPrePostProcessors/`.
- Supports both explicit and implicit stepping paths.
- Tracks `activationTime` when `Vm` crosses threshold.
- Exports selected ionic variables to volumetric fields.

Stimulus architecture:

- PDE stimulus is loaded from `monodomainStimulus` dictionary entries.
- Ionic-model S1/S2 stimulus can exist for single-cell workflows.
- Guard logic prevents accidental double stimulation in monodomain runs unless explicitly allowed (`allowIonicStimulusInMonodomain`).

### 3) `SingleCellSolver`

Single integration-point workflow:

- Instantiates ionic model with one integration point.
- Enables `solveVmWithinODESolver=true` so voltage is advanced inside ODEs.
- Writes traces to `postProcessing/<ionicModel>_<tissue>_<stimulus>.txt`.
- Uses shared `ionicModelIO` helpers for headers and row writes.

### 4) `EikonalSolver`

Reduced-order activation model:

- Solves steady eikonal-diffusion formulation for activation time `psi`.
- Uses anisotropic conductivity tensor and optional stabilized form (`eikonalAdvectionDiffusionApproach`).
- Applies stimulus through geometric box selection in the mesh.

### 5) Optional staged domains

`electrophysicsSystem` can optionally own staged auxiliary domains:

- A pre-primary conduction-system domain, such as Purkinje, advanced before
  the myocardium and coupled through a domain-coupling model.
- A post-primary ECG domain, advanced after the myocardium from the updated
  electrophysiology state with one-way data flow.

Current ECG implementations live under `ecgModels/`. The
`PseudoECGDomain` model remains a one-way ECG evaluation, but it is now
scheduled through the same domain system as the Purkinje side.

Canonical dictionary schema for staged domains (inside
`MonoDomainSolverCoeffs` / `BiDomainSolverCoeffs`):

```cpp
conductionNetworkDomains
{
    purkinjeNetwork
    {
        ConductionSystemDomain purkinjeNetworkModel;
        // ... purkinjeNetworkModel settings ...
    }
}

domainCouplings
{
    purkinjeToMyocardium
    {
        ElectroDomainCoupler pvjResistanceCouplingModel;
        // ... coupling settings ...
    }
}

ecgDomains
{
    ECG
    {
        ECGDomain PseudoECGDomain;
        // ... ECG model settings ...
    }
}
```

Only the schema above is supported.

### 6) `ElectroMechanicalModel`

Wrapper workflow for coupled electro-mechanics runs:

- Owns an electro model plus a runtime-selected active-tension model.
- Passes coupling signals from the ionic model stack into the active-tension
  stack.
- Uses the same runtime-selection pattern as the other electro models.

## Runtime chain with `cardiacFoam`

1. `cardiacFoam` creates `physicsModel`.
2. Electro `physicsModel` creates `electroModel` from `electroProperties`.
3. `MonoDomainSolver` and `SingleCellSolver` create `ionicModel` from the
   selected electro-model coeff dictionary in `electroProperties`.

## Build target

`Make/files` builds:
all compiled sources under:

- `core/electroModel/`
- `domains/`
- `conductionSystemModels/`
- `solvers/`
- `wrappers/`

into:

- `$(FOAM_USER_LIBBIN)/libelectroModels`

Electro-domain coupling models are compiled from `src/couplingModels/electroDomain/`.
