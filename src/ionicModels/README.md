# ionicModels

The cell models: published cardiac ionic models that give the electrical response of a cell. The myocardium and Purkinje domains each run one, and single-cell runs use them directly. Every model is converted from its CellML description by [cellML2foam](../../applications/scripts/cellML2foam/README.md).

## What's available

Grouped by the tissue each model represents. "Biophysical" models describe the individual ionic currents; "phenomenological" models reproduce the action potential with a few variables.

### Ventricular

| Model | Species | Type | Reference |
|---|---|---|---|
| `TNNP` | human | biophysical | ten Tusscher, Noble, Noble and Panfilov (2004), American Journal of Physiology-Heart and Circulatory Physiology, [doi:10.1152/ajpheart.00794.2003](https://doi.org/10.1152/ajpheart.00794.2003) |
| `BuenoOrovio` | human | phenomenological, four variables | Bueno-Orovio, Cherry and Fenton (2008), Journal of Theoretical Biology, [doi:10.1016/j.jtbi.2008.03.029](https://doi.org/10.1016/j.jtbi.2008.03.029) |
| `ToRORd_dynCl` | human | biophysical, with dynamic intracellular chloride | Tomek et al. (2020), bioRxiv, [doi:10.1101/2020.06.01.127043](https://doi.org/10.1101/2020.06.01.127043) |
| `TWorld` | human | biophysical, with calcium handling, contraction and β-adrenergic signalling | Tomek et al. (2025), bioRxiv, [doi:10.1101/2025.03.24.645031](https://doi.org/10.1101/2025.03.24.645031) |
| `Gaur` | pig | biophysical | Gaur et al. (2021), PLOS Computational Biology, [doi:10.1371/journal.pcbi.1009137](https://doi.org/10.1371/journal.pcbi.1009137) |
| `AlievPanfilov` | dog | phenomenological, two variables | Aliev and Panfilov (1996), Chaos, Solitons & Fractals, [doi:10.1016/0960-0779(95)00089-5](https://doi.org/10.1016/0960-0779(95)00089-5) |

### Atrial

| Model | Species | Type | Reference |
|---|---|---|---|
| `Courtemanche` | human | biophysical | Courtemanche, Ramirez and Nattel (1998), American Journal of Physiology-Heart and Circulatory Physiology, [doi:10.1152/ajpheart.1998.275.1.H301](https://doi.org/10.1152/ajpheart.1998.275.1.H301) |
| `Grandi` | human | biophysical, sinus rhythm and chronic atrial fibrillation | Grandi et al. (2011), Circulation Research, [doi:10.1161/CIRCRESAHA.111.253955](https://doi.org/10.1161/CIRCRESAHA.111.253955) |
| `PerisYague` | pig | biophysical, extends Courtemanche with chloride currents | Peris-Yagüe et al. (2022), Frontiers in Physiology, [doi:10.3389/fphys.2022.812535](https://doi.org/10.3389/fphys.2022.812535) |

### Sinoatrial node

| Model | Species | Type | Reference |
|---|---|---|---|
| `Fabbri` | human | biophysical, pacemaking | Fabbri et al. (2017), The Journal of Physiology, [doi:10.1113/JP273259](https://doi.org/10.1113/JP273259) |

### Purkinje

| Model | Species | Type | Reference |
|---|---|---|---|
| `Stewart` | human | biophysical | Stewart et al. (2009), Philosophical Transactions of the Royal Society A, [doi:10.1098/rsta.2008.0283](https://doi.org/10.1098/rsta.2008.0283) |
| `Trovato` | human | biophysical, automaticity | Trovato et al. (2020), Journal of Molecular and Cellular Cardiology, [doi:10.1016/j.yjmcc.2020.04.001](https://doi.org/10.1016/j.yjmcc.2020.04.001) |

### Batched versions

Every model also comes as a GPU-ready batched version: `AlievPanfilovcompactBatched`, `BuenoOroviocompactBatched`, `CourtemanchecompactBatched`, `FabbricompactBatched`, `GaurcompactBatched`, `GrandicompactBatched`, `PerisYaguecompactBatched`, `StewartcompactBatched`, `TNNPcompactBatched`, `TWorldcompactBatched`, `ToRORd_dynClcompactBatched` and `TrovatocompactBatched`. Batched models run on the CPU; building with `CARDIAC_ENABLE_CUDA` set adds a CUDA path.

### Cell types

`TNNP`, `BuenoOrovio`, `ToRORd_dynCl` and `TWorld` carry their own endocardial, mid-myocardial and epicardial constants. The other eight accept those labels as well, but their equations hold one cell type, so a label there only picks which `ionicConstantOverrides` scope applies: without overrides every region gets identical cells. In batched form only the same four support the three cell types; the other eight run `myocyte` cells only. Cell types can vary across the tissue, by transmural bands, named regions or mesh cell zones, with an optional gradient on top.

### Verification models

Three manufactured-solution models, used by [verificationModels](../verificationModels/README.md): `monodomainFDAManufactured`, `bidomainFDAManufactured` and `bathBidomainFDAManufactured`.

## Folders

```text
src/ionicModels/
├── ionicModel/          # base classes, batched base, heterogeneity support
├── verificationModels/  # the three manufactured-solution models
├── <Name>/              # one folder per model, with its generated equations
└── <Name>Batched/       # one folder per batched model
```

**Deep dive:** [IONIC_MODEL_ARCHITECTURE.md](IONIC_MODEL_ARCHITECTURE.md) explains how this library is built inside.

## What this does not own

- The domains that run these models: [electroModels](../electroModels/README.md).
- Active tension: [activeTensionModels](../activeTensionModels/README.md), which reads `Vm` or `Cai` through `ElectromechanicalSignalProvider`, implemented by `ionicModel`.
- The verifiers: [verificationModels](../verificationModels/README.md).
