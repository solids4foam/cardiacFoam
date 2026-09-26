# NiedererEtAl2011 tutorial architecture

This tutorial implements the Niederer slab verification workflow for tissue-scale
monodomain simulations.

- Electro model: `monodomainSolver`
- Typical ionic model: `TNNP` (corrected 2026-09-26: this said `BuenoOrovio`; the committed `constant/electroProperties` sets `ionicModel TNNP;`)
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
    ionicModel TNNP;    // corrected 2026-09-26, matches constant/electroProperties
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

Driver-managed (from the repository root; a record stages into a scratch
directory you supply, never into this tree):

```bash
driverFoam run --strict --entry niederer2011 --cases-root tutorials --scratch-dir <dir>
driverFoam sweep-run --spec tutorials/NiedererEtAl2011verification/setup/studies/cartesianConvergence/sweep_hex_convergence.json --output-dir <dir>
driverFoam sweep-run --spec tutorials/NiedererEtAl2011verification/setup/studies/tetConvergence/sweep_tet_generic.json --output-dir <dir>
```

The driver reads this case as it is: `niederer2011` is a pointer at this
directory, and each study under `setup/studies/` states only what it
varies (`dx`/`tetDx`, `deltaT`, `endTime`). Each study's `base` names
`cases_root` (`tutorials`, relative to the repository root), because a
driver sweep over a tutorial record has no cases root it could discover.

## Regression behavior

`regression/regressionTest.sh` is local to this case and validates against
`regression/NiedererEtAl2011.reference`.
