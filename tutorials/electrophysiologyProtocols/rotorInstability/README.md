# rotorInstability tutorial architecture

This tutorial implements a 2D slab workflow demonstrating self-sustained rotor
dynamics. It uses multiple stimuli to initiate and sustain spiral waves (rotors)
in the tissue.

- Electro model: `monodomainSolver`
- Typical ionic model: `BuenoOrovio`
- Main metric: activation-time behavior and self-sustained reentry

## Folder structure

```text
tutorials/coreProtocols/rotorInstability/
├── constant/
│   ├── electroProperties
│   └── polyMesh/ (generated)
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   ├── fvSolution
│   └── rotorpoints
├── rotorInstability.reference
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
2. runs `cardiacFoam` (serial or parallel)
3. runs OpenFOAM postProcess function objects for probe extraction using `rotorpoints`

## Run modes

Manual:

```bash
./Allrun
./Allrun parallel
./regressionTest.sh
```

## Regression behavior

`regressionTest.sh` is local to this case and validates against the reference
file `rotorInstability.reference` stored in the case directory.
