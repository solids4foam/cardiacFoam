# cardiacFoam

`cardiacFoam` is an OpenFOAM toolbox for cardiac electrophysiology and electro-mechanics.
It keeps the `solids4foam` runtime-selection style, and can run in two modes:

- Full mode: with a full `solids4foam` installation (electro + solid + FSI workflows).
- Electro-only mode: with the lightweight fallback `modules/physicsModel` shipped in this repository.

## Repository architecture

```text
cardiacFoam/
├── applications/
│   ├── solvers/cardiacFoam/                 # Main executable
│   ├── utilities/                           # sweepCurrents, setFibreField, ...
│   └── scripts/
│       ├── driverFoam/openfoam_driver/      # Python tutorial automation engine
│       └── cellML2foam/                     # CellML -> ionic model generation pipeline
├── src/
│   ├── electroModels/                       # Runtime-selectable electro solvers
│   ├── ionicModels/                         # Runtime-selectable ionic ODE models
│   ├── genericWriter/                       # Shared I/O and stimulus parsing helpers
│   └── couplingModels/                      # Electro-mechanics coupling signal interfaces
├── modules/
│   └── physicsModel/                        # Lightweight fallback physicsModel
├── tutorials/                               # Reference cases and regression cases
└── etc/resolveSolids4Foam.sh                # Backend selection helper
```

## Runtime execution flow

At runtime, solver/model selection is fully dictionary-driven:

1. `applications/solvers/cardiacFoam/cardiacFoam.C` creates `physicsModel::New(runTime)`.
2. `physicsModel` type is selected from `constant/physicsProperties` (`type`).
3. For electro runs, `src/electroModels/electroModel::New(...)` selects `electroModel` from `constant/electroProperties` (`electroModel`).
4. Electro models (`MonoDomainSolver`, `SingleCellSolver`, `EikonalSolver`) select ionic models through `ionicModel::New(...)` (`ionicModel` in electro coefficients).

This gives one stable executable (`cardiacFoam`) with pluggable electro and ionic sub-models.

## Current electro model stack

- `MonoDomainSolver`: tissue PDE-ODE model, explicit/implicit stepping, activation-time tracking.
- `SingleCellSolver`: single integration-point ODE workflow (no spatial PDE solve).
- `EikonalSolver`: reduced-order activation-time model.

## Current ionic model stack

Compiled in `libionicModels`:

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
- `tmanufacturedFDA` (manufactured-solution verification model)

## Tutorial and automation architecture

Manual tutorial entry points are in `tutorials/*/Allrun`. For parameter sweeps and post-processing automation, use:

- `applications/scripts/driverFoam/openfoam_driver`

Current tutorial specs in the Python driver:

- `singleCell`
- `niederer2012`
- `manufacturedFDA`
- `restitutionCurves`

The driver writes run manifests and artifact manifests (`run_manifest.json`, `plots.json`) for reproducibility.

## Build and run

```bash
./Allwmake
```

Run examples:

```bash
cd tutorials/singleCell
./Allrun

cd ../NiedererEtAl2012
./Allrun parallel
```

Run automation from repository root:

```bash
foamctl all --tutorial singleCell
foamctl all --tutorial niederer2012
```

## Regression checks

- Tutorial regression entrypoint: `tutorials/Alltest-regression`
- Driver architecture contract tests: `applications/scripts/driverFoam/openfoam_driver/tests/test_tutorial_architecture_contract.py`

## Notes

This is active research software. APIs and dictionaries evolve as models and workflows are expanded.
