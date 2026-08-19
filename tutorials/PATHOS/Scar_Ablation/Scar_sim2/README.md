# Scar_sim2

Monodomain tissue simulation with a scar region modelled through reduced
conductivity in the `conductivity` field. The case uses a Purkinje network
coupled to the myocardium and computes pseudo-ECG output.

## Stack

- electro model: `monodomainSolver`
- ionic model: `BuenoOrovio`
- scar: reduced `conductivity` field in the scar cellZone (set in mesh pre-processing)
- conduction system: `purkinjeGraphModel` with `monodomain1DSolver`
- ECG: `ecgDomains.ECG.ecgSolver pseudoECG`

## Dictionary scope

`constant/electroProperties`:

```cpp
myocardiumSolver monodomainSolver;

monodomainSolverCoeffs
{
    ionicModel  BuenoOrovio;
    tissue      epicardialCells;

    ionicHeterogeneity { ... }

    solutionAlgorithm implicit;

    conductionNetworkDomains { purkinjeNetwork { ... } }
    domainCouplings           { purkinjeToMyocardium { ... } }

    ecgDomains
    {
        ECG { ecgSolver pseudoECG; electrodePositions { ... } }
    }
}
```

## Outputs

- `postProcessing/pseudoECG.dat`

## Execution

Mesh must be pre-converted (VTK import + `transformPoints`) before running:

```bash
./Allrun
./Allrun parallel
```
