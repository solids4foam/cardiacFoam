# electroDomains

The domains the core builds. Each one owns the state of one physical region and hands the numerics to a solver.

## What's available

| Domain | What it holds |
|---|---|
| `myocardiumDomain` | the heart tissue: `Vm`, its ionic model and a tissue solver (monodomain or bidomain). `eikonalMyocardiumDomain` is its eikonal form. |
| `conductionSystemDomain` | the Purkinje network as a graph, with its ionic model and a Purkinje solver |
| `ecgDomain` | electrodes and an ECG model; it reads the tissue state and never changes it |
| `extracellularPotentialDomain` | the global extracellular potential `phiE` over heart and bath, shared with a bidomain myocardium |

The domains exchange data only through interfaces: coupling endpoints, used by the Purkinje–muscle junction couplers, and state providers, used by the ECG.

## Folders

```text
src/electroModels/electroDomains/
├── myocardiumDomain/
├── conductionSystemDomain/
├── ecgDomain/
└── extracellularPotentialDomain/
```

The solvers these domains use are in [myocardiumModels](../myocardiumModels/README.md), [conductionSystemModels](../conductionSystemModels/README.md) and [ecgModels](../ecgModels/README.md); the couplers are in [electroCouplers](../electroCouplers/README.md).
