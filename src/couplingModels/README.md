# couplingModels

The link between the two kinds of cell-level model. It holds headers only; there is no compiled library.

## What's available

- `electromechanicalSignalProvider.H`: the `ElectromechanicalSignalProvider` contract. `ionicModel` implements it and `activeTensionModel` reads from it, so a tension model gets `Vm` or `Cai` without depending on a particular cell model.

## Folders

```text
src/couplingModels/
└── electromechanicalSignalProvider.H
```

## What this does not own

- The couplers between electrical domains (Purkinje, ECG, bath): [electroModels/electroCouplers](../electroModels/electroCouplers/README.md).
