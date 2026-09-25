# restitutionCurves_s1s2Protocol tutorial architecture

This tutorial runs single-cell S1-S2 pacing sweeps to generate restitution curves.

- Electro model: `singleCellSolver`
- Typical ionic model: configurable (default driver uses `BuenoOrovio`)
- Purpose: APD/restitution analysis across S2 intervals

## Folder structure

```text
tutorials/coreProtocols/restitutionCurves_s1s2Protocol/
├── constant/
│   ├── electroProperties
│   ├── physicsProperties
│   └── sweepCurrents
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── setup/
│   ├── run_cases.sh
│   ├── setup_multiple_simulations_s1s2.py
│   ├── postProcessing_restCurves.py
│   ├── animate_trace.py
│   └── mainRestitutionCurves_s1s2Protocol.py
├── plotVoltage
├── Allrun
├── Allclean
└── README.md
```

## Key dictionary scope

`constant/electroProperties`:

```cpp
myocardiumSolver singleCellSolver;

singleCellSolverCoeffs
{
    ionicModel ...;
    tissue ...;

    singleCellStimulus
    {
        stim_period_S1 ...;
        nstim1 ...;
        stim_period_S2 ...;
        nstim2 ...;
    }
}
```

## Run modes

Manual:

```bash
./Allrun
```

omniD-managed sweeps (this case IS the default -- omniD's `restitutionCurves`
tutorial record points at this directory and stages a clone of it for every
case; nothing here is ever written in place):

```bash
omnidriver --plugin cardiacfoam describe --entry restitutionCurves --cases-root <path to tutorials>
omnidriver --plugin cardiacfoam sweep-plan --spec setup/studies/tworldS1S2Restitution/sweep.json --output-dir <scratch output dir>
omnidriver --plugin cardiacfoam sweep-run  --spec setup/studies/tworldS1S2Restitution/sweep.json --output-dir <scratch output dir>
```

`setup/studies/tworldS1S2Restitution/sweep.json` is this tutorial's own
study, in omniD's tutorial-record vocabulary (docs/superpowers/specs/
2026-09-24-tutorials-are-pointers-design.md in the omniD repository): a
study name is either a literal `document:dotted.path` dictionary key
(`constant/electroProperties:singleCellSolverCoeffs.tissue`, set directly)
or one of this record's two allowed axes -- `ionicModel` (a bare model
name; derives `singleCellSolverCoeffs.ionicModel` and this model's catalogued
single-cell `stim_amplitude`) and `s1s2Protocol` (a mapping of
`s1_interval_ms`/`n_s1`/`s2_interval_ms`/`n_s2`; derives the
`singleCellStimulus` S1/S2 keys plus the case's `endTime`/`writeAfterTime`).
This study reproduces the sweep the old, now-deleted `driver_config.json`
named (TWorld, S1=1000ms/10 beats, S2 swept from 1500ms down to 250ms) --
`driver_config.json` itself is gone: it spoke to the retired driverFOAM
add-on's own Python defaults, which no longer exist.
