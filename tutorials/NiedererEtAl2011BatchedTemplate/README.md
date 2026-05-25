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

## CPU vs GPU (Batched) Comparisons

This tutorial includes scripts to benchmark the adaptive CPU ODE solver (`RKF45`) against the fixed-step batched integrators (`batched_euler`, `batched_heun`, `batched_rl`, `batched_soa`).

To run the full suite of comparisons:

```bash
./runComparisons.sh
```

This master script performs two steps sequentially:
1. **Mode Comparison:** Compares all models at a baseline configuration (`run_niederer_bueno_orovio_batched_comparison.sh`).
2. **Substep Sweep:** Runs `euler` and `soa` at different substep counts (5, 10, 20, 50, 100) to trace out the speed-vs-accuracy tradeoff (`run_substep_sweep.sh`).

### Where to find results

All outputs are generated in the `comparisonResults/` folder:

- **Metrics Tables**: `comparison_metrics.json` and `substep_sweep_metrics.json`
- **Plots (`comparisonResults/plots/`)**:
  - `bar_point_errors.png`: Per-point absolute errors across models.
  - `line_activation_profile.png`: Arc-length vs activation time along a diagonal line.
  - `line_abs_error_vs_arclength.png`: Absolute error of activation time along the diagonal line.
  - `scatter_points_wall_time_vs_max_error.png`: Speed vs accuracy scatter plot.
  - `substep_*.png`: Sweeps demonstrating convergence and tradeoff as the number of substeps changes.
