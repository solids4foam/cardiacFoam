# BrugadaSyndrome

Monodomain tissue simulation modelling Brugada syndrome type-1 pattern via
`ionicConstantOverrides` on the TNNP ionic model, with pseudo-ECG output and
a coupled Purkinje network.

## Stack

- electro model: `monodomainSolver`
- ionic model: `TNNP`
- ionic overrides: INa −80% (global), Ito ×6 + ICaL −50% (epicardialCells only)
- conduction system: `purkinjeGraphModel` with `monodomain1DSolver`
- ECG: `ecgDomains.ECG.ecgSolver pseudoECG`

## Dictionary scope

`constant/electroProperties`:

```cpp
myocardiumSolver monodomainSolver;

monodomainSolverCoeffs
{
    ionicModel  TNNP;
    tissue      epicardialCells;

    ionicHeterogeneity { ... }

    ionicConstantOverrides
    {
        global          { scale { g_Na 0.2; } }
        epicardialCells { scale { g_to 6.0; g_CaL 0.5; } }
    }

    conductionNetworkDomains { purkinjeNetwork { ... } }
    domainCouplings           { purkinjeToMyocardium { ... } }

    ecgDomains
    {
        ECG
        {
            ecgSolver pseudoECG;
            electrodePositions { V1 ...; V2 ...; ... }
        }
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
