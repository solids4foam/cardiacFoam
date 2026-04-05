# couplingModels architecture

This directory now groups shared coupling contracts used across electro-only and
electro-mechanics workflows.

## Current contents

```text
src/couplingModels/
├── common/                     # Cross-library coupling interfaces
│   └── electromechanicalSignalProvider.H
├── electroDomain/              # Electro-domain coupling contracts + models
│   ├── electroDomainCouplingEndpoints.H
│   ├── ElectroDomainCoupler.{H,C}
│   └── pvj*.{H,C}
├── lnInclude/
└── README.md
```

## `CouplingSignalProvider` interface (`common/`)

Purpose:

- Provide a minimal API for querying cell-level coupling signals from
  electrophysiology models.

Main methods:

```cpp
virtual bool hasSignal(CouplingSignal s) const = 0;
virtual scalar signal(label i, CouplingSignal s) const = 0;
```

`Foam::ionicModel` inherits from this interface; active-tension models consume
it via `electroMechanicalModel`.

## Electro-domain coupling contracts (`electroDomain/`)

`electroDomainCouplingEndpoints.H` and `ElectroDomainCoupler.H` define
the coupling interface between staged electrophysiology domains (e.g. Purkinje
to myocardium). Runtime models such as `pvjResistanceCouplingModel` implement
that interface.
