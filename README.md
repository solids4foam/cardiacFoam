# cardiacFoam

`cardiacFoam` is an OpenFOAM toolbox for cardiac electrophysiology and electromechanics. It implements PDE-ODE solvers for tissue-level propagation, single-cell ODE integration, eikonal activation models, ECG computation, and operator-split electromechanical coupling.

It runs in two modes:

- **Full mode** — with a `solids4foam` installation (`modules/solids4foam`). Enables electromechanical coupling and solid mechanics workflows.
- **Electro-only mode** — with the lightweight `modules/physicsModel` shipped in this repository. No solid mechanics dependency required.

## Repository layout

```text
cardiacFoam/
├── src/
│   ├── electroModels/          # Electro solvers (monodomain, bidomain, eikonal, ECG, conduction system)
│   ├── ionicModels/            # Ionic ODE models (serial and GPU-batched variants)
│   ├── activeTensionModels/    # Active tension models for electromechanical coupling
│   ├── electroMechanicalModels/# Operator-split electromechanical coupling (solids4foam builds only)
│   ├── couplingModels/         # Signal interfaces between electro and solid solvers
│   ├── verificationModels/     # Manufactured-solution verification utilities
│   └── genericWriter/          # Shared I/O and stimulus parsing
├── applications/
│   ├── solvers/cardiacFoam/    # Main solver executable
│   ├── utilities/              # Pre/post-processing utilities (mesh, fibres, ECG, Purkinje, ...)
│   └── scripts/
│       ├── driverFoam/         # Python automation engine for tutorials and parameter sweeps
│       └── cellML2foam/        # CellML → ionic model code generation pipeline
├── modules/
│   ├── physicsModel/           # Lightweight physicsModel fallback (electro-only builds)
│   └── solids4foam/            # solids4foam submodule (full builds)
├── tutorials/                  # Reference and research cases
└── etc/resolveSolids4Foam.sh   # Build-mode selection helper
```

Each subdirectory carries its own `README.md` (and `ARCHITECTURE.md` where relevant) with component-level detail.

## Runtime architecture

For a spatial electrophysiology case, the selection and ownership flow is:

```text
cardiacFoam
  -> physicsModel::New()                 constant/physicsProperties
  -> electroModel::New()                 constant/electroProperties
  -> electrophysiologyModel              multi-domain orchestration
  -> myocardiumDomain
  -> myocardiumSolver                    PDE kernel
```

`physicsModel` selects `electroModel` through the `type` entry.
`electroModel` then reads the public `myocardiumSolver` key. The spatial runtime
names `monodomainSolver`, `bidomainSolver`, and `eikonalSolver` all enter the
`electrophysiologyModel` assembly path; `singleCellSolver` is a direct
`electroModel` implementation and bypasses the multi-domain builder. Domain
objects own fields and state, solver objects own numerical kernels, and `core/`
owns orchestration. See [`src/electroModels/ARCHITECTURE.md`](src/electroModels/ARCHITECTURE.md).

## What the code contains

**Electro solvers** (`src/electroModels/`) — dictionary-driven runtime selection across four domains: myocardium PDE-ODE (monodomain, bidomain), eikonal activation, ECG forward problem, and 1D conduction-system models (Purkinje, restitution-aware eikonal). See [`src/electroModels/README.md`](src/electroModels/README.md).

**Ionic models** (`src/ionicModels/`) — a library of human and animal cardiac cell models. Every model ships a serial variant and a GPU-batched variant for tissue-scale simulations. Manufactured-solution verification models are included for FDA-style solver validation. See [`src/ionicModels/README.md`](src/ionicModels/README.md).

**Active tension models** (`src/activeTensionModels/`) — active stress generation models (Nash–Panfilov, Goktepe–Kuhl, Land–Niederer) with serial and GPU-batched variants, and a manufactured-solution verification layer for coupled electromechanics. See [`src/activeTensionModels/README.md`](src/activeTensionModels/README.md).

**Electromechanical coupling** (`src/electroMechanicalModels/`, `src/couplingModels/`) — sequential operator-split coupling of the electro and solid solvers. Requires a solids4foam build. See [`src/electroMechanicalModels/README.md`](src/electroMechanicalModels/README.md).

**Utilities** (`applications/utilities/`) — mesh and fibre setup, Purkinje graph runner, ECG recomputation, ionic heterogeneity probing, current sweep, VTK conversion, and more. Each utility has its own README.

**Driver** (`applications/scripts/driverFoam/`) — a Python automation engine for running tutorials, parameter sweeps, and post-processing. Produces reproducible run and artifact manifests. See [`applications/scripts/driverFoam/openfoam_driver/README.md`](applications/scripts/driverFoam/openfoam_driver/README.md).

**Tutorials** (`tutorials/`) — organised into three groups:

| Group | Contents |
|---|---|
| `electrophysiologyProtocols/` | Single-cell ODE runs, restitution curves, heterogeneity probes, rotor dynamics |
| `manufacturedSolutions/` | MMS verification cases for all solver variants incl. electromechanics |
| `NiedererEtAl2011/` | Benchmark cases: tissue propagation, Purkinje, electromechanics |

See [`tutorials/README.md`](tutorials/README.md).

## Build

```bash
./Allwmake
```

Requires an OpenFOAM environment. For GPU-batched ionic models, see [`src/ionicModels/IONIC_MODEL_ARCHITECTURE.md`](src/ionicModels/IONIC_MODEL_ARCHITECTURE.md). For build-mode selection (full vs electro-only), see [`etc/resolveSolids4Foam.sh`](etc/resolveSolids4Foam.sh).

Electromechanical tutorials require a full solids4foam build and will not run in electro-only mode.

The repository owns the code under `src/`, `applications/`, the lightweight
`modules/physicsModel` fallback, and the CellML generator/templates. Generated
ionic-model equation headers should be changed through that generation path.
`modules/solids4foam` is an external submodule and is not maintained as
cardiacFoam source.

## Regression

```bash
tutorials/Alltest-regression
```

## Notes

This is active research software. APIs, dictionaries, and model interfaces evolve as the toolbox is extended.
