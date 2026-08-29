# TWORLD versus Gaur species comparison

This study compares the pig Gaur model with the human TWORLD model in the
single-cell solver using the Land--Niederer active-tension model.

## Contents

- `sweep_tworld_vs_gaur.json`: six-case zip sweep across model/species and
  pacing cycle length;
- `postprocess_tworld_vs_gaur.py`: postprocessor for waveforms, rate dependence, and
  metrics;
- `results/<run-name>/`: one generated run directory containing the manifest,
  `cases/`, figures, metrics, and optional summaries.

The intended layout is:

```text
tworldVsGaur/
├── sweep_tworld_vs_gaur.json
├── postprocess_tworld_vs_gaur.py
├── README.md
└── results/
    └── <run-name>/
        ├── sweep_manifest.json
        ├── cases/
        ├── species_comparison_waveforms.png
        ├── species_comparison_all_variables.png
        ├── species_comparison_rate_dependence.png
        └── species_comparison_metrics.csv
```

## Run

From the repository root, source the host OpenFOAM installation and run:

```bash
source /Volumes/OpenFOAM-v2412/etc/bashrc
export DRIVERFOAM_RUNTIME_CONFIG=/Users/simaocastro/omnidriver/driverfoam-runtime.yaml

applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/sweep_tworld_vs_gaur.json \
    --output-dir tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/results/sweepRun

python3 tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/postprocess_tworld_vs_gaur.py \
    --input-dir tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/results/sweepRun/cases \
    --output-dir tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/results/sweepRun
```

The study covers CL = 300, 500, and 1000 ms. The focused waveform PNG is a
2×1 presentation figure: Vm with faded Ca²⁺ on the right-hand y-axis above,
then Ta alone. The detailed PNG and raw traces also
include ICaL, IKr, IK1, Ito, model-specific SR release, and SERCA fluxes.
The transparent `species_comparison_calcium_overlay.png` contains only the
faded Ca²⁺ curves and their right-hand axis for overlaying on another voltage
figure. Each run directory contains the driver manifest, raw case outputs, the
PNG figures, and the metrics CSV. Use a new descriptive run name for each
fresh execution; all generated files under `results/` are ignored by git.
