# activeTensionModels library architecture

This directory provides runtime-selectable active-tension models built into
`libactiveTensionModels`.

These models sit on the electromechanics side of the stack: they consume
coupling signals from electrophysiology models and evolve active-tension state
variables for each integration point.

## Directory structure

```text
src/activeTensionModels/
├── activeTensionModel/   # Base class and runtime selection
├── verificationModels/
│   └── ManufacturedElectromechanics/ # Manufactured electromechanics verification model
├── GoktepeKuhl/          # Goktepe-Kuhl phenomenological active tension model
├── NashPanfilov/         # Nash-Panfilov phenomenological active tension model
├── LandNiederer/         # Original Land et al. intact-human model
├── LandNiedererTWorld/   # TWorld-derived six-state contraction subsystem
├── *Batched/             # GPU-ready Batched versions of the models (e.g. NashPanfilovBatched)
├── Make/
├── lnInclude/
└── README.md
```

## Core class: `Foam::activeTensionModel`

Defined in `activeTensionModel/activeTensionModel.H` and implemented in
`activeTensionModel/activeTensionModel.C`.

Main responsibilities:

- Runtime selection via `activeTensionModel::New(...)`.
- Store the active-tension dictionary and integration-point count.
- Query an optional `ElectromechanicalSignalProvider` for upstream signals such as `Vm`
  and `Cai`.
- Run the base `calculateTension(...)` loop and delegate per-point work to
  `solveAtPoint(...)`.
- Provide shared I/O and export helpers through `activeTensionIO`.

## Batched Execution (`Foam::batchedActiveTensionModel`)

Models extending `batchedActiveTensionModel` use a Structure-of-Arrays (SoA) data layout
to dispatch execution across CPU (OpenMP) or GPU (CUDA) backends. The base
`calculateTension(...)` loop is overridden to target the selected backend,
mirroring the `ionicModels` batched pattern.

## Runtime model contract

`libactiveTensionModels` is built in both modes, but spatial electromechanical
workflows that consume these models require full solids4foam mode. The values in
the table are the exact registered `activeTensionModel` dictionary selectors.

| Purpose | Runtime name | Backend/build availability | Source boundary | Known equivalence limitation |
|---|---|---|---|---|
| Phenomenological tension | `GoktepeKuhl` | scalar CPU; full EM workflows | maintained wrapper; generated equations/Names metadata | Batched integration uses a different data path |
| Phenomenological tension | `GoktepeKuhlBatched` | SoA host, optional CUDA; full EM workflows | maintained wrapper/backend; generated batch equations | Scalar/batched trajectories require tolerance-based comparison |
| Phenomenological tension | `NashPanfilov` | scalar CPU; full EM workflows | maintained wrapper; generated equations/Names metadata | Batched integration uses a different data path |
| Phenomenological tension | `NashPanfilovBatched` | SoA host, optional CUDA; full EM workflows | maintained wrapper/backend; generated batch equations | Scalar/batched trajectories require tolerance-based comparison |
| Biophysical tension | `LandNiederer` | scalar CPU; full EM workflows | original seven-state intact-human model | Active output is `AV_Ta`; passive and total tension remain diagnostic outputs |
| Biophysical tension | `LandNiedererTWorld` | scalar CPU; full EM workflows | TWorld six-state contraction subsystem | Applies resting-Cai preconditioning |
| Biophysical tension | `LandNiedererTWorldBatched` | SoA host, optional CUDA; full EM workflows | TWorld batched backend | Uses explicit batched resting-Cai conditioning; scalar equivalence requires stated tolerances |
| Electromechanical MMS | `ManufacturedElectromechanics` | scalar CPU; full EM verification | maintained verification implementation and Names metadata | Verification-only; not a physiological tension law |

All production models select their driving electrophysiology signal from
dictionary input (`couplingSignal`, normally `Vm` or `Cai`) and integrate with
the `ElectromechanicalSignalProvider` interface used by `ionicModel`.
Both scalar Land variants use `preconditioningTime` (default 1000 ms) and
integrate to the resting steady state over a fixed 100 substeps.
`LandNiedererTWorldBatched` advances its generated hot path, conditions one
representative state, and copies that state to every integration point.

### Generated and maintained boundaries

- Runtime wrappers, registration, coupling-signal handling, backend dispatch,
  preconditioning, and CUDA kernels listed by `Make/files-gpu` are maintained.
- Model equation and Names headers are generated inputs. Update their generator
  contract and regenerate instead of editing equations directly.
- CUDA adds a backend to the three registered batched selectors; it does not add
  new selector names. Host and CUDA results require explicit parity tolerances.
- `NiedererHunterSmith` is commented out in `Make/files`; it is not compiled or
  runtime-selectable and is therefore not advertised as available.
