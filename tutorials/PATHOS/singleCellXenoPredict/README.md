# singleCell tutorial architecture

This case is the single integration-point electrophysiology workflow.

- Electro model: `singleCellSolver`
- Voltage evolution: inside ionic ODE system (`solveVmWithinODESolver=true`)
- Spatial PDE solve: not used

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

## Run modes

Manual:

```bash
./Allrun
./regressionTest.sh
```

Driver-managed sweep:

```bash
foamctl all --entry singleCell
```

The Python driver mutates ionic model, tissue, and stimulus amplitude for each case,
then collects outputs and post-processes in `setup`.

## Regression behavior

`regressionTest.sh` is local to this case and validates against
`singleCell.reference`.
