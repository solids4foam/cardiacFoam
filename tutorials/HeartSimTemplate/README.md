# NiedererEtAl2012 tutorial architecture

This tutorial implements the Niederer slab verification workflow for tissue-scale
monodomain simulations.

- Electro model: `MonoDomainSolver`
- Typical ionic model: `TNNP`
- Main metric: activation-time behavior and smoke-check fields

## Folder structure

```text
tutorials/NiedererEtAl2012/
├── 0/
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
├── Allrun
├── runRegressionTest.sh
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

`runRegressionTest.sh` executes the smoke-check function object and compares against
`system/smokeCheck.reference`.
