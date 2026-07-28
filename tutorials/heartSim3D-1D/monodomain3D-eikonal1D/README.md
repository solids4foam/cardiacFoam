# AcuteIschemia

Monodomain tissue simulation modelling early acute ischemia via
`ionicConstantOverrides` on the TNNP ionic model, with pseudo-ECG output and
a coupled Purkinje network.

## Stack

- electro model: `monodomainSolver`
- ionic model: `TNNP`
- ionic overrides: hyperkalemia (K_o 5.4 → 8.0 mM), INa −40%, ICaL −40%, IK1 ×2, INaK −50%
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
        global
        {
            set   { K_o 8.0; }
            scale { g_Na 0.6; g_CaL 0.6; g_K1 2.0; P_NaK 0.5; }
        }
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
