# singleCell tutorial architecture

This case is the single integration-point electrophysiology workflow.

- Electro model: `singleCellSolver`
- Voltage evolution: inside ionic ODE system (`solveVmWithinODESolver=true`)
- Spatial PDE solve: not used

## Folder structure

```text
tutorials/coreProtocols/singleCell/
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
│   ├── run_cases.sh
│   └── singleCellinteractivePlots.py
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
applications/scripts/driverFoam/bin/driverFoam sweep-plan \
    --spec tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/sweep_tworld_vs_gaur.json \
    --output-dir tutorials/electrophysiologyProtocols/singleCell/setup/studies/tworldVsGaur/results/sweepRun
applications/scripts/driverFoam/bin/driverFoam sweep-run \
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

Driver-managed sweep:

```bash
applications/scripts/driverFoam/bin/driverFoam run --strict --entry singleCell
```

The Python driver mutates ionic model, tissue, and stimulus amplitude for each case,
then collects outputs and post-processes in `setup`.

## Regression behavior

`regressionTest.sh` is local to this case and validates against
`singleCell.reference`.
