# NiedererEtAl2011 tutorial architecture

This tutorial implements the Niederer slab verification workflow for tissue-scale
monodomain simulations.

- Electro model: `monodomainSolver`
- Typical ionic model: `BuenoOrovio`
- Main metric: activation-time behavior and smoke-check fields

## Folder structure

```text
tutorials/NiedererEtAl2011/
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   ├── fvSolution
│   └── smokeCheck.reference
├── setupNiedererEtAl2011/
│   └── postProcessing/
│       ├── cache_postProcessing.py
│       ├── line_postProcessing.py
│       └── points_postProcessing.py
├── NiedererEtAl2012.reference
├── regressionTest.sh
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
./regressionTest.sh
./regressionTest.sh parallel
```

Driver-managed sweep:

```bash
foamctl all --entry niederer2012
```

Driver sweeps are controlled by
`applications/scripts/driverFoam/openfoam_driver/core/defaults/niederer_2012.py`.

## Regression behavior

`regressionTest.sh` is local to this case and validates against the reference
file stored in the case directory.
