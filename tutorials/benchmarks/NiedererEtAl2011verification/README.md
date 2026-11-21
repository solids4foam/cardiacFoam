# NiedererEtAl2011 tutorial architecture

## What it does

This tutorial implements the Niederer et al. (2011) slab benchmark for
tissue-scale monodomain simulations: a 20 x 3 x 7 mm slab, the TNNP ionic model
(`TNNPcompactBatched` with `batchedIntegrator rushLarsen`) and the benchmark's
corner stimulus. It writes `activationTime` at the nine `system/Niedererpoints`
probes and along the `system/Niedererlines` diagonal.

- **Exercises:** the monodomain solver driving a real activation front from a
  published ionic model through 3D tissue, which the smooth manufactured cases
  do not have; the batched TNNP model with Rush-Larsen gating.
- **Use it for:** reproducing the Niederer et al. (2011) activation-time
  benchmark and its space/time convergence grid (`setup/studies/`), and as a
  starting point for slab studies with a published ionic model.

## Regression

Tier: `full`.

`regression/regressionTest.sh` runs `./Allrun parallel` once (6 subdomains,
`scotch`) with the case's own settings: 100 x 15 x 35 hexahedra (0.2 mm),
`deltaT 1e-5`, `endTime 0.015`. It compares six `activationTime` probes of
`Niedererpoints` at 15 ms against `regression/NiedererEtAl2011.reference`
(6 rows, tolerance 1e-4 s): the stimulated corner (1.19 ms) and the
`(0, 3, 7) mm` corner (11.2 ms) have activated, and the far end
`x = 20 mm` (three corners) and the slab centre still read -1. The values
are frozen from an earlier run of the code, not the published benchmark
values.

```bash
./regression/regressionTest.sh
```

## Folder structure

```text
tutorials/benchmarks/NiedererEtAl2011verification/
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
│   ├── line_postProcessing.py
│   ├── points_postProcessing.py
│   ├── table_summary.py
│   └── studies/
│       ├── cartesianConvergence/
│       └── tetConvergence/
├── regression/
│   ├── regressionTest.sh
│   ├── NiedererEtAl2011.reference
│   └── full-only
├── bidomain/          the same benchmark with bidomainSolver
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
    ionicModel TNNPcompactBatched;
    batchedIntegrator rushLarsen;
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
3. runs `postProcess` with the `Niedererpoints` and `Niedererlines` probe
   function objects on the latest time

## Run modes

Manual:

```bash
./Allrun
./Allrun parallel
./regression/regressionTest.sh
```

Omnidriver-managed (from the repository root; a record stages into a scratch
directory you supply, never into this tree):

```bash
omnidriver --plugin cardiacfoam --repo . plan --strict \
    --entry niederer2011 --scratch-dir <scratch>
```

omnidriver reads this case as it is: `niederer2011` is a pointer at this
directory, and each study under `setup/studies/` states only what it
varies (`dx`/`tetDx`, `deltaT`, `endTime`). Each study's `base` names
`cases_root` (`tutorials`, relative to the repository root), because a
omnidriver sweep over a tutorial record has no cases root it could discover.

`cartesianConvergence` and `tetConvergence` run Niederer et al. (2011)'s
grid: Δx = 0.5, 0.2, 0.1 mm (`dx`/`tetDx` 5e-4, 2e-4, 1e-4 m) × Δt = 0.05,
0.01, 0.005 ms (`deltaT` 5e-5, 1e-5, 5e-6 s), nine cases each, with `endTime`
0.2, 0.08 and 0.055 s for Δx 0.5, 0.2 and 0.1 mm, chosen so every probe has
activated.
The Δx 0.1 mm cases take hours each and are meant for an HPC run.

## Notes on dictionary entries

`constant/electroProperties` selects the batched TNNP model with Rush-Larsen
gating: the same wavefront as TNNP, much faster. The batched models do not use
the `ODESolver`, so there is no `solver` or `maxSteps` entry.
