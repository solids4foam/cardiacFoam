# restitutionCurves_s1s2Protocol tutorial architecture

## What it does

Runs single-cell S1-S2 pacing (`singleCellSolver`, one integration point) to
generate action-potential restitution curves. The committed case paces a
Bueno-Orovio epicardial cell ten times at a 2000 ms cycle length, follows
with two S2 beats at 250 ms, integrates with the Euler ODE solver, and writes
the Vm trace after 16 s to `postProcessing/`; the
`tworldS1S2Restitution` study sweeps the S2 interval for TWorld.

- **Exercises:** `singleCellSolver` with the S2 pulses of
  `singleCellStimulus` (`stim_period_S2`, `nstim2`), which no other case
  sets, and the Euler ODE integrator.
- **Use it for:** APD and restitution curves of an ionic model across S2
  intervals (`setup/postProcessing_restCurves.py`), before the model is used
  in tissue.

## Regression

Tier: `none`.

None. `singleCellSolver` and its trace output are checked by
[`singleCell`](../singleCell/README.md) (standard), and the Bueno-Orovio
ionic model by [`ionicHeterogeneity`](../ionicHeterogeneity/README.md)
(standard), which integrates it in the `ionicHeterogeneityProbe` utility.
What this case adds is the restitution measurement protocol. The S2 branch
of `singleCellStimulus` and the Euler ODE integrator are not run by any
regression.

## Folder structure

```text
tutorials/electrophysiologyProtocols/restitutionCurves_s1s2Protocol/
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
│   ├── postProcessing_restCurves.py
│   ├── table_summary.py
│   ├── animate_trace.py
│   └── studies/
│       └── tworldS1S2Restitution/
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

omniD-managed sweeps (this case IS the default -- omniD's `restitutionCurves`
tutorial record points at this directory and stages a clone of it for every
case; nothing here is ever written in place):

```bash
omnidriver --plugin cardiacfoam describe --entry restitutionCurves --cases-root <path to tutorials>
omnidriver --plugin cardiacfoam sweep-plan --spec setup/studies/tworldS1S2Restitution/sweep.json --output-dir <scratch output dir>
omnidriver --plugin cardiacfoam sweep-run  --spec setup/studies/tworldS1S2Restitution/sweep.json --output-dir <scratch output dir>
```

`setup/studies/tworldS1S2Restitution/sweep.json` is this tutorial's own
study, in omniD's tutorial-record vocabulary: a
study name is either a literal `document:dotted.path` dictionary key
(`constant/electroProperties:singleCellSolverCoeffs.tissue`, set directly)
or one of this record's two allowed axes -- `ionicModel` (a bare model
name; derives `singleCellSolverCoeffs.ionicModel` and this model's catalogued
single-cell `stim_amplitude`) and `s1s2Protocol` (a mapping of
`s1_interval_ms`/`n_s1`/`s2_interval_ms`/`n_s2`; derives the
`singleCellStimulus` S1/S2 keys plus the case's `endTime`/`writeAfterTime`).
This study sweeps TWorld at S1=1000ms/10 beats with S2 from 1500ms down to 250ms.

## Mesh resolutions

`system/blockMeshDict` spans 20 x 3 x 7 mm with (200 30 70) cells, 0.1 mm.
`singleCellSolver` integrates one ionic cell whatever the mesh, so the cell
count does not enter the solution.
