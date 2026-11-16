# restitutionCurves_s1s2Protocol tutorial architecture

## What it does

Runs single-cell S1-S2 pacing (`singleCellSolver`, one integration point) to
generate action-potential restitution curves. The committed case paces a
Bueno-Orovio epicardial cell ten times at a 2000 ms cycle length, follows
with two S2 stimuli at 250 ms, integrates with the Euler ODE solver, and writes
the Vm trace after 16 s to `postProcessing/`; the
`tworldS1S2Restitution` study sweeps the S2 interval for TWorld.

- **Exercises:** `singleCellSolver` with the S2 pulses of
  `singleCellStimulus` (`stim_period_S2`, `nstim2`), which no other case
  sets, and the Euler ODE integrator, which no other regression uses.
- **Use it for:** APD and restitution curves of an ionic model across S2
  intervals (`setup/postProcessing_restCurves.py`), before the model is used
  in tissue.

## Regression

Tier: `standard`.

`regression/regressionTest.sh` runs the tutorial as it is (serial, about
10 s) and checks `Vm` of
`postProcessing/BuenoOrovio_epicardialCells_S1_2000_S2_250.txt` at four
times (`regression/restitutionCurves_s1s2Protocol.reference`, 4 rows,
tolerance 5e-3 mV, values frozen from an earlier run). The first S2 stimulus,
250 ms after the last S1 beat, falls in its refractory tail and does not
capture; the second does. The rows sample the tail of the last S1 beat
(18.25 s), the blocked S2 (18.3 s), and the plateau and repolarisation of the
captured S2 beat (18.6 s, 18.7 s). Run it alone with:

```bash
./regression/regressionTest.sh
```

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
├── regression/
│   ├── regressionTest.sh
│   └── restitutionCurves_s1s2Protocol.reference
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

## Mesh

`system/blockMeshDict` is one hexahedron, as in `../singleCell`:
`singleCellSolver` integrates one ionic cell whatever the mesh.
