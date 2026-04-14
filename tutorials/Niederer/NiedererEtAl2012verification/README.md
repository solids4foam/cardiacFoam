# NiedererEtAl2012 tutorial architecture

This tutorial implements the Niederer slab verification workflow for tissue-scale
monodomain simulations.

- Electro model: `MonoDomainSolver`
- Typical ionic model: `TNNP`
- Main metric: activation-time behavior and smoke-check fields

## Folder structure

```text
tutorials/NiedererEtAl2012/
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
├── setupNiedererEtAl2012/
│   └── postProcessing/
│       ├── cache_postProcessing.py
│       ├── line_postProcessing.py
│       └── points_postProcessing.py
├── runRegressionTest.sh
├── Allrun
├── Allclean
└── README.md
```

## Key dictionary scope

`constant/electroProperties`:

```cpp
electroModel MonoDomainSolver;

MonoDomainSolverCoeffs
{
    ionicModel TNNP;
    tissue epicardialCells;
    solutionAlgorithm implicit;   // or explicit via sweeps

    monodomainStimulus
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
./runRegressionTest.sh
./runRegressionTest.sh parallel
```

Driver-managed sweep:

```bash
foamctl all --tutorial niederer2012
```

Driver sweeps are controlled by
`applications/scripts/driverFoam/openfoam_driver/core/defaults/niederer_2012.py`.

## Regression behavior

`runRegressionTest.sh` is a local wrapper around
`tutorials/regressionTests/NiedererEtAl2012/runRegressionTest.sh`. The shared
regression folder stores the override dictionaries and
`NiedererEtAl2012.reference`, while the runnable case stays focused on the base
tutorial setup.
