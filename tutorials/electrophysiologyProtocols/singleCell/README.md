# singleCell tutorial architecture

This case is the single integration-point electrophysiology workflow.

- Electro model: `singleCellSolver`
- Voltage evolution: inside ionic ODE system (`solveVmWithinODESolver=true`)
- Spatial representation: single integration point (0-D)

## Folder structure

```text
tutorials/electrophysiologyProtocols/singleCell/
├── constant/
│   ├── electroProperties
│   ├── physicsProperties
│   └── sweepCurrents
├── system/
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── setup/
│   ├── singleCellinteractivePlots.py
│   ├── table_summary.py
│   ├── sweep_ionic_model_tissue.json
│   └── studies/
│       └── tworldVsGaur/
├── singleCell.reference
├── regressionTest.sh
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
        stim_start ...;
        stim_duration ...;
        stim_amplitude ...;
        stim_period_S1 ...;
        nstim1 ...;
        stim_period_S2 ...;
        nstim2 ...;
    }
}
```

## Outputs

`singleCellSolver` writes traces to:

- `postProcessing/<ionicModel>_<tissue>_<stimulusSuffix>.txt`

Optional plotting is done by `plotVoltage` (skipped by default when `CF_SKIP_PLOTS=1`).

## TWORLD versus Gaur species comparison

The study is kept under `setup/studies/tworldVsGaur/`, following the repository
study convention:

```text
setup/studies/tworldVsGaur/
├── sweep_tworld_vs_gaur.json
├── postprocess_tworld_vs_gaur.py
└── results/                 # generated and gitignored
```

`sweep_tworld_vs_gaur.json` defines the focused comparison used for
pig versus human ventricular cells:

- Gaur / `myocyte` as the pig case;
- TWORLD / `endocardialCells` as the human case;
- pacing cycle lengths of 1000, 500, and 300 ms;
- `activeTensionModel LandNiedererTWorld` for every case.

The exported traces include Vm, `cai`, ICaL, Jrel/Jup (using each model's
native names), IKr, IK1, Ito, and the TWorld contraction `AV_Ta` trace. Run the
study from the repository root with:

```bash
driverFoam sweep-plan \
    --spec tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/sweep_tworld_vs_gaur.json \
    --output-dir tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/results/sweepRun
driverFoam sweep-run \
    --spec tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/sweep_tworld_vs_gaur.json \
    --output-dir tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/results/sweepRun
python3 tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/postprocess_tworld_vs_gaur.py \
    --input-dir tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/results/sweepRun/cases \
    --output-dir tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/results/sweepRun
```

This writes the focused 2×1 `species_comparison_waveforms.png`, the detailed
`species_comparison_all_variables.png`,
`species_comparison_rate_dependence.png`, and
`species_comparison_calcium_overlay.png`, and
`species_comparison_metrics.csv` into the run directory. The waveform figure
uses the 1000-ms beat and contains Vm alone above, followed by Ta alone; the
separate transparent calcium overlay carries the Ca²⁺ curves and right-hand
y-axis. The detailed figure also
shows ICaL, SR fluxes, and repolarisation currents. The transparent calcium
overlay contains only the Ca²⁺ traces and their right-hand axis. The rate figure reports
APD90 and peak Ca versus pacing cycle length. Keep each completed run in a
separate `results/<run-name>/` directory.

## Run modes

Manual:

```bash
./Allrun
./regressionTest.sh
```

omniD-managed sweeps (this case IS the default -- omniD's `singleCell`
tutorial record points at this directory and stages a clone of it for every
case; nothing here is ever written in place):

```bash
omnidriver --plugin cardiacfoam describe --entry singleCell --cases-root <path to tutorials>
omnidriver --plugin cardiacfoam sweep-plan --spec setup/sweep_ionic_model_tissue.json --output-dir <scratch output dir>
omnidriver --plugin cardiacfoam sweep-run  --spec setup/sweep_ionic_model_tissue.json --output-dir <scratch output dir>
```

`setup/sweep_ionic_model_tissue.json` and `setup/studies/tworldVsGaur
/sweep_tworld_vs_gaur.json` are this tutorial's own studies, in omniD's
tutorial-record vocabulary (docs/superpowers/specs/2026-09-24-tutorials-
are-pointers-design.md in the omniD repository): a study name is either a
literal `document:dotted.path` dictionary key
(`constant/electroProperties:singleCellSolverCoeffs.tissue`, set directly)
or this record's one allowed axis, `ionicModel` (a bare model name; derives
`singleCellSolverCoeffs.ionicModel` and this model's catalogued single-cell
`stim_amplitude`). Neither study varies S2 pacing, so there is no
`s1s2Protocol`-style axis here (`restitutionCurves`'s own axis of that name
derives `endTime`/`writeAfterTime` from all four S1+S2 numbers together --
nothing this case's studies vary needs deriving); `stim_period_S1` and
`outputVariables.ionic.export` are set as direct keys too.

`setup/driver_config.json`, which named a single ionic model for the now-
retired driverFOAM add-on's own Python defaults, is deleted: it named no
key any current path reads. `sweep_tworld_vs_gaur.json`'s old `base` also
carried `stim_start_ms`/`n_s1`/`end_time_buffer_s`/`write_after_time_s`,
none of which the old Python factory's `make_spec` ever accepted as
keywords -- dead values with no effect before this rewrite either, now
dropped rather than carried forward as new dead keys.

## Regression behavior

`regressionTest.sh` is local to this case and validates against
`singleCell.reference`.
