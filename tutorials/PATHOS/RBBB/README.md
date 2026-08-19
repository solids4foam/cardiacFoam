# RBBB

Monodomain tissue simulation modelling right bundle branch block (RBBB) via a
modified Purkinje graph with the right bundle severed. No `ionicConstantOverrides`
are used — the conduction delay is entirely structural.

## Stack

- electro model: `monodomainSolver`
- ionic model: `BuenoOrovio`
- conduction system: `purkinjeGraphModel` with `monodomain1DSolver`, graph `purkinjeGraph.rbbb`
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

    conductionNetworkDomains
    {
        purkinjeNetwork
        {
            conductionSystemDomain  purkinjeGraphModel;
            purkinjeGraphModelCoeffs
            {
                conductionSystemSolver monodomain1DSolver;
                graphFile              purkinjeGraph;   // symlinked to purkinjeGraph.rbbb
                ...
            }
        }
    }

    domainCouplings { purkinjeToMyocardium { ... } }

    ecgDomains
    {
        ECG { ecgSolver pseudoECG; electrodePositions { ... } }
    }
}
```

## Outputs

- `postProcessing/pseudoECG.dat`

## Execution

```bash
./Allrun
./Allrun parallel
```
