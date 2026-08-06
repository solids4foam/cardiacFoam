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

Driver-managed sweeps:

```bash
foamctl all --entry restitutionCurves
```

Driver defaults live in
`applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/defaults/restitution_curves.py`.

The driver mutates ionic model/tissue/stimulus values per case, updates end time,
collects `.txt` outputs, and can generate per-case animations before post-processing.
