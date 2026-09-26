# cardiacFoam

[![Build and test](https://github.com/solids4foam/cardiacFoam/actions/workflows/buildAndTest.yml/badge.svg)](https://github.com/solids4foam/cardiacFoam/actions/workflows/buildAndTest.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

cardiacFoam is an OpenFOAM toolbox for cardiac electrophysiology and electromechanics. It covers the range from a single cell to the ventricles: action potentials from 12 published cell models, activation across tissue (monodomain, bidomain or eikonal), Purkinje networks, ECGs, and, with solids4foam, active tension and contraction. Every model is chosen by name in the case's dictionaries.

It combines two parts: the electrophysiology libraries in this repository and [solids4foam](https://github.com/solids4foam/solids4foam) for the solid mechanics. Electrophysiology runs on its own; electromechanics is electrophysiology plus solids4foam.

## Get it

```bash
git clone https://github.com/solids4foam/cardiacFoam.git
cd cardiacFoam
```

You need OpenFOAM **v2312, v2406, v2412, v2506 or v2512**, sourced in your shell. Electromechanics also needs solids4foam; see [Build modes](#build-modes).

## Build

```bash
./Allwmake
```

## Run a first case

```bash
cd tutorials/electrophysiologyProtocols/singleCell
./Allrun
```

This runs a single TWorld endocardial cell, paced twice at 1000 ms, for 2 s of simulated time. The voltage trace is written to `postProcessing/TWorld_endocardialCells_S1_1000.txt`, and `./Allclean` resets the case. All the other cases are listed in [tutorials](tutorials/README.md).

## Simulate the ventricles

`tutorials/idealizedHeart/electroHeart` simulates electrophysiology on an idealized biventricular geometry (123,617 cells), with a Purkinje network and a pseudo-ECG. Its mesh is stored with Git LFS, and the case reads `../mesh`, so copy the whole `idealizedHeart` folder and set the run length on your copy:

```bash
git lfs install && git lfs pull
mkdir -p "$FOAM_RUN" && cp -r tutorials/idealizedHeart "$FOAM_RUN"/
cd "$FOAM_RUN"/idealizedHeart/electroHeart
foamDictionary system/controlDict.monodomain -entry endTime -set 0.6
./Allrun parallel
```

The tutorial stops at 0.04 s so that its regression test stays short; the `foamDictionary` line sets a longer run, 0.6 s here. At a `deltaT` of 2e-5 that is 30,000 steps, so `./Allrun parallel` splits it over 6 processors, which takes roughly 10 minutes on a recent laptop. Leave out `parallel` to run in serial.

## Where to go next

- [tutorials/](tutorials/README.md): runnable cases, verification studies and benchmarks.
- [src/](src/README.md): the libraries and what each one is for.
- [applications/utilities/](applications/utilities/README.md): tools for preparing and post-processing cases.
- [src/verificationModels/](src/verificationModels/README.md) and [tutorials/manufacturedSolutions/](tutorials/manufacturedSolutions/README.md): how the numerics are verified.

## Build modes

`./Allwmake` chooses the mode through `etc/resolveSolids4Foam.sh`:

- **Full.** Used when `SOLIDS4FOAM_INST_DIR` points at a built solids4foam or, if that variable is unset, when a built copy is found in `modules/solids4foam`, `~/solids4foam` or `$WM_PROJECT_USER_DIR/solids4foam`. Everything builds, including electromechanics.
- **Lightweight.** Used when no built solids4foam is found, or when you set `FORCE_LIGHTWEIGHT_PHYSICSMODEL=1`. Everything except electromechanics builds.

To use the bundled solids4foam:

```bash
git submodule update --init --recursive
(cd modules/solids4foam && ./Allwmake)
```

## Testing

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight ./tutorials/Alltest-regression
```

For a full build, use `CARDIAC_REGRESSION_BUILD_MODE=with-solids4foam`. The mode must match how you built.

To run a subset, pass substrings of the case paths, e.g.
`./tutorials/Alltest-regression bidomain purkinje`; `-l` lists the cases.
Each case's `regression/regressionTest.sh` can use the shared helpers in
`tutorials/regressionFunctions`, e.g. `checkSolverLogs` to fail when a solver
run did not finish cleanly.

## Licence and citation

cardiacFoam is active research software, released under the [GNU GPL v3](LICENSE). If you use it in published work, please cite: **cardiacFOAM**.
