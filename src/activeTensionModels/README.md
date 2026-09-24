# activeTensionModels

Cell-level models of the mechanical response: the same kind of model as the ionic ones, turning voltage or calcium into an active tension `Ta` for the solid. The electromechanics part of the core runs them, and single-cell runs can compute tension too.

## What's available

| Model | What it is | Reference |
|---|---|---|
| `NashPanfilov` / `NashPanfilovBatched` | phenomenological tension, driven by `Vm` | Nash and Panfilov (2004), Progress in Biophysics and Molecular Biology, [doi:10.1016/j.pbiomolbio.2004.01.016](https://doi.org/10.1016/j.pbiomolbio.2004.01.016) |
| `LandNiederer` / `LandNiedererBatched` | the seven-state intact-human model, driven by `Cai` | Land et al. (2017), Journal of Molecular and Cellular Cardiology, [doi:10.1016/j.yjmcc.2017.03.008](https://doi.org/10.1016/j.yjmcc.2017.03.008) |
| `LandNiedererTWorld` / `LandNiedererTWorldBatched` | the six-state contraction subsystem of T-World, driven by `Cai` | Tomek et al. (2025), bioRxiv, [doi:10.1101/2025.03.24.645031](https://doi.org/10.1101/2025.03.24.645031) |
| `ManufacturedElectromechanics` | manufactured-solution model for verification, driven by `Vm` | — |

Each model reads its one signal through `ElectromechanicalSignalProvider` in [couplingModels](../couplingModels/README.md). The batched versions run on the CPU, with an optional CUDA path.

## Folders

```text
src/activeTensionModels/
├── activeTensionModel/   # base classes and shared helpers
├── verificationModels/   # ManufacturedElectromechanics
├── NashPanfilov/
├── LandNiederer/
├── LandNiedererTWorld/
└── <Name>Batched/        # the batched version of each model
```

**Deep dive:** [ARCHITECTURE.md](ARCHITECTURE.md) explains how this library is built inside.

## What this does not own

- The passive mechanical response: the solid owns it. `LandNiederer`'s passive branch is diagnostic output only.
- The signals it reads: [couplingModels](../couplingModels/README.md).
- When tension is computed: [electroMechanicalModels](../electroMechanicalModels/README.md).
