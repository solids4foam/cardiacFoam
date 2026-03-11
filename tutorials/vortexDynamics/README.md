# vortexDynamics tutorial architecture

This tutorial is a tissue-scale 2D propagation/reentry-style setup using the
monodomain workflow.

- Electro model: `monoDomainElectro`
- Ionic model: configured in `constant/electroProperties`
- Purpose: wave dynamics demonstration in a planar domain

## Folder structure

```text
tutorials/vortexDynamics/
├── 0/
│   └── Vm
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── Allrun
├── Allclean
└── README.md
```

## Key dictionary scope

`constant/electroProperties` uses `monoDomainElectroCoeffs` with:

- conductivity, `chi`, `Cm`
- `ionicModel` and `tissue`
- `monodomainStimulus` (supports multiple stimulus boxes and start times)
- solution settings (`solutionAlgorithm`, ODE solver controls)

## Run modes

```bash
./Allrun
./Allrun parallel
```

`Allrun` builds mesh (`blockMesh`) and then runs `cardiacFoam` in serial or parallel.
