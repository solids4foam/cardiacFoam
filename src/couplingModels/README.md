# couplingModels

This folder contains shared coupling-side interfaces used across libraries.
It does not contain the staged Purkinje, ECG, or bath electro-domain couplers.

## Current contents

```text
src/couplingModels/
├── common/
│   └── electromechanicalSignalProvider.H
├── lnInclude/
└── README.md
```

## Purpose

The main role of this folder is to define the
`ElectromechanicalSignalProvider` contract used by:

- `ionicModel`
- `activeTensionModel`

That interface lets one component query scalar signals such as:

- `Vm`
- `Cai`

without depending on a concrete ionic-model implementation.

## What this folder does not own

This folder does not own:

- Purkinje-to-myocardium electro couplers
- ECG or bath couplers
- staged domain-coupling runtime selection

Those pieces live in:

- [`../electroModels/electroCouplers/README.md`](../electroModels/electroCouplers/README.md)

