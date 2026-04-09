# manufacturedSolutions/singleCellBidomain tutorial

This case is the single integration-point manufactured-solution workflow for the
bidomain ionic stack.

- Electro model: `SingleCellSolver`
- Voltage evolution: inside ionic ODE system (`solveVmWithinODESolver=true`)
- Spatial PDE solve: not used

## Folder structure

```text
tutorials/manufacturedSolutions/singleCellBidomain/
├── constant/
│   ├── electroProperties
│   ├── physicsProperties
│   └── sweepCurrents
├── system/
│   ├── blockMeshDict
│   └── controlDict
├── setupSingleCell/
│   ├── run_cases.sh
│   └── singleCellinteractivePlots.py
├── Allrun
├── runRegressionTest.sh
├── Allclean
└── README.md
```

## Key dictionary scope

`constant/electroProperties`:

```cpp
electroModel SingleCellSolver;

SingleCellSolverCoeffs
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

`SingleCellSolver` writes traces to:

- `postProcessing/<ionicModel>_<tissue>_<stimulusSuffix>.txt`

Optional plotting is done by `plotVoltage` (skipped by default when `CF_SKIP_PLOTS=1`).

## Run modes

Manual:

```bash
./Allrun
./runRegressionTest.sh
```

Driver-managed sweep:

```bash
applications/scripts/driverFoam/bin/driverFoam all --tutorial singleCell --config tutorials/manufacturedSolutions/singleCellBidomain/setupSingleCell/driver_config.json
```

The Python driver mutates ionic model, tissue, and stimulus amplitude for each case,
then collects outputs and post-processes in `setupSingleCell`.

## Regression behavior

`runRegressionTest.sh` validates selected points from output traces against
`system/singleCell.reference`.
