# cardiacFoam

`cardiacFoam` is an OpenFOAM toolbox for cardiac electrophysiology and electro-mechanics.
It can run in two modes:

- Full mode: with a full `solids4foam` installation (electro + solid + FSI workflows).
- Electro-only mode: with the lightweight fallback `modules/physicsModel` shipped in this repository.

## Repository architecture

```text
cardiacFoam/
├── applications/
│   ├── solvers/cardiacFoam/                 # Main executable
│   ├── utilities/                           # sweepCurrents, setFibreField, runPurkinjeGraph, ...
│   └── scripts/
│       ├── driverFoam/openfoam_driver/      # Python tutorial automation engine
│       └── cellML2foam/                     # CellML → ionic model generation pipeline
├── src/
│   ├── electroModels/                       # Runtime-selectable electro solvers
│   ├── ionicModels/                         # Runtime-selectable ionic ODE models
│   ├── activeTensionModels/                 # Runtime-selectable active-tension models
│   ├── verificationModels/                  # Manufactured-solution verification layer
│   ├── electroMechanicalModels/             # Electro-mechanics coupling (solids4foam builds only)
│   ├── genericWriter/                       # Shared I/O and stimulus parsing helpers
│   └── couplingModels/                      # Electro-mechanics coupling signal interfaces
├── modules/
│   └── physicsModel/                        # Lightweight fallback physicsModel
├── tutorials/                               # Reference cases and regression cases
└── etc/resolveSolids4Foam.sh                # Backend selection helper
```

## Runtime execution flow

At runtime, solver/model selection is fully dictionary-driven:

1. `applications/solvers/cardiacFoam/cardiacFoam.C` creates the run time.
2. `applications/solvers/cardiacFoam/cardiacFoam.C` creates `physicsModel::New(runTime)`.
3. `physicsModel` type is selected from `constant/physicsProperties` (`type`).
4. For electro runs, `src/electroModels/electroModel::New(...)` selects the solver from `constant/electroProperties` (`myocardiumSolver`).
5. Electro models (`monodomainSolver`, `bidomainSolver`, `singleCellSolver`, `eikonalSolver`) select ionic models through `ionicModel::New(...)` (`ionicModel` in electro coefficients).

## Current electro model stack

- `monodomainSolver`: tissue PDE-ODE model, explicit/implicit stepping, activation-time tracking.
- `bidomainSolver`: coupled intra- and extracellular tissue PDE-ODE model.
- `singleCellSolver`: single integration-point ODE workflow (no spatial PDE solve).
- `eikonalSolver`: reduced-order activation-time model.

## Current ionic model stack

Compiled in `libionicModels`:

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
- `monodomainFDAManufactured` (manufactured monodomain verification model)
- `bidomainFDAManufactured` (manufactured bidomain verification model)
- `bathBidomainFDAManufactured` (manufactured bath-bidomain verification model)

## Tutorial and automation architecture

Manual tutorial entry points are in `tutorials/*/Allrun`. For parameter sweeps and post-processing automation, use:

- `applications/scripts/driverFoam/openfoam_driver`

Current tutorial specs in the Python driver:

- `singleCell`
- `niederer2012`
- `manufacturedFDA`
- `manufacturedFDABidomain`
- `manufacturedFDABathBidomain`
- `manufacturedEikonalECG`
- `manufacturedMonodomainTotalLagrangianEM`
- `restitutionCurves`

The driver writes run manifests and artifact manifests (`run_manifest.json`, `plots.json`) for reproducibility.

## Build and run

```bash
./Allwmake
```

Run examples:

```bash
cd tutorials/singleCellprotocols/singleCell
./Allrun

cd tutorials/NiedererEtAl2011/NiedererEtAl2011verification
./Allrun parallel
```

Run automation from repository root:

```bash
foamctl all --entry singleCell
foamctl all --entry niederer2012
```

Electromechanical tutorials such as
`tutorials/NiedererEtAl2011/electroMechanicalNiedererEtAl2011` require a full
`solids4foam`-backed build of `cardiacFoam`. They do not run in the lightweight
electro-only mode.

## Regression checks

- Tutorial regression entrypoint: `tutorials/Alltest-regression`
- Driver architecture contract tests: `applications/scripts/driverFoam/openfoam_driver/tests/test_tutorial_architecture_contract.py`

## Notes

This is active research software. APIs and dictionaries evolve as models and workflows are expanded.
