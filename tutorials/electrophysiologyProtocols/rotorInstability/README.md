# rotorInstability tutorial architecture

## What it does

Solves the monodomain equation (`monodomainSolver`, explicit) with
Bueno-Orovio epicardial cells on a 50 x 50 x 0.1 mm slab and uses three
stimuli to start and sustain a spiral wave (rotor): an S1 strip on the left
edge at 0 s, an S2 on the lower half at 0.45 s, and a strip on the right edge
at 2 s. It runs to 4 s and samples `activationTime` at five probes
(`system/rotorpoints`).

- **Exercises:** `monodomainSolver` with the explicit algorithm and a
  three-entry `externalStimulus` list, `activationTime` updated by repeated
  re-activation, and the probe function object. Only this case checks
  sustained re-entry: activations after the last stimulus that come from the
  rotor alone.
- **Use it for:** a starting point for re-entry and rotor studies on a 2D
  slab; checking that an ionic model or conductivity change still sustains
  (or terminates) a rotor.

## Regression

Tier: `full` (marker `regression/full-only`).

`regression/regressionTest.sh` runs `Allrun parallel` (6 subdomains,
`scotch`) with `endTime` 1 instead of the tutorial's 4 s, so the third
stimulus at 2 s is not reached. It compares the last activation time at each
of the five probes at t = 1.0 s (5 rows in
`regression/rotorInstability.reference`) with frozen values from an earlier
run, within 2e-3 s. By 1 s every probe has been activated again by the rotor
that the S2 stimulus starts at 0.45 s, so the check requires sustained
re-entry. Run it on its own from the case root:

```bash
./regression/regressionTest.sh
```

## Folder structure

```text
tutorials/electrophysiologyProtocols/rotorInstability/
├── constant/
│   ├── electroProperties
│   ├── physicsProperties
│   └── polyMesh/ (generated)
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   ├── fvSolution
│   └── rotorpoints
├── regression/
│   ├── full-only
│   ├── regressionTest.sh
│   └── rotorInstability.reference
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
    solutionAlgorithm explicit;

    externalStimulus
    {
        // Three distinct stimuli to induce reentry:
        stimulusLocationMinList ((0 0 0) (0 0 0) (40e-3 0 0));
        stimulusLocationMaxList ((2e-3 50e-3 0.1e-3) (50e-3 25e-3 0.1e-3) (50e-3 50e-3 0.1e-3));
        stimulusDurationList (2e-3 5e-3 2e-3);
        stimulusIntensityList (50000 50000 50000);
        stimulusStartTimeList (0.0 0.45 2.0);
    }
}
```

## Execution flow

`Allrun`:

1. runs `blockMesh`
2. runs `cardiacFoam` (serial, or in parallel with `decomposePar`,
   `reconstructPar` and `redistributeRestartState`)
3. runs OpenFOAM postProcess function objects for probe extraction using `rotorpoints`

## Run modes

Manual:

```bash
./Allrun
./Allrun parallel
./regression/regressionTest.sh
```

## Mesh

`system/blockMeshDict` spans 50 x 50 x 0.1 mm with 150 x 150 x 1 cells,
deltaX = 0.33 mm.
