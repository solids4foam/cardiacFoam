# NiedererEtAl2011 tutorial architecture

This tutorial implements the Niederer slab verification workflow for tissue-scale
monodomain simulations.

- Electro model: `monodomainSolver`
- Typical ionic model: `BuenoOrovio`
- Main metric: activation-time behavior and smoke-check fields

## Folder structure

```text
tutorials/NiedererEtAl2011verification/
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   ├── fvSolution
│   ├── Niedererlines
│   └── Niedererpoints
├── setup/
│   ├── convert_raw_samples.py
│   ├── line_postProcessing.py
│   ├── points_postProcessing.py
│   ├── table_summary.py
│   └── studies/
│       ├── cartesianConvergence/
│       └── tetConvergence/
├── regression/
│   ├── regressionTest.sh
│   └── NiedererEtAl2011.reference
├── Allrun
├── Allclean
└── README.md
```

## Key dictionary scope

`constant/electroProperties`:

```cpp
myocardiumSolver monodomainSolver;

monodomainSolverCoeffs
{
    ionicModel BuenoOrovio;
    tissue epicardialCells;
    solutionAlgorithm implicit;   // or explicit via sweeps

    externalStimulus
    {
        ...
    }
}
```

## Execution flow

`Allrun`:

1. runs `blockMesh`
2. runs `cardiacFoam` (serial or parallel)
3. runs OpenFOAM postProcess function objects for smoke checks and probe extraction

## Run modes

Manual:

```bash
./Allrun
./Allrun parallel
regression/regressionTest.sh
regression/regressionTest.sh parallel
```

Driver-managed sweep:

```bash
driverFoam run --strict --entry niederer2012
```

Driver sweeps are controlled by the driverFOAM add-on's `niederer_2012`
cardiacFoam plugin defaults.

## Regression behavior

`regression/regressionTest.sh` is local to this case and validates against
`regression/NiedererEtAl2011.reference`.
