# electroCouplers

This directory contains inter-domain coupling infrastructure for the staged
electrophysiology system. It defines the endpoint contracts used by domains,
the base class for couplers, and the Purkinje-ventricular-junction (PVJ)
coupling family.

## Current contents

```text
src/electroModels/electroCouplers/
├── electroDomainCouplingEndpoints.H   # Domain-side coupling interfaces
├── electroDomainCoupler.{H,C}         # Base class for staged couplers
├── heartBathInterfaceCoupler.{H,C}    # Bath-interface coupling code still present
├── pvjCoupler/
│   ├── pvjCoupler.{H,C}               # PVJ family base class
│   ├── pvjMapper.{H,C}                # PVJ geometry and source projection
│   ├── eikonal/
│   │   └── eikonalPvjCoupler.{H,C}
│   └── reactionDiffusion/
│       └── reactionDiffusionPvjCoupler.{H,C}
└── README.md
```

## Endpoint contracts

`electroDomainCouplingEndpoints.H` defines the minimal interfaces shared by
coupling models:

- `tissueCouplingEndpoint`
  - implemented by `MyocardiumDomain`
  - exposes the tissue mesh, `Vm`, and the volumetric `sourceField`
- `networkCouplingEndpoint`
  - implemented by `ConductionSystemDomain`
  - exposes terminal-node locations and voltages
  - accepts terminal coupling currents and matching diagnostic source terms

These interfaces keep the domain classes decoupled from any specific coupling
implementation.

## Base coupling model

`ElectroDomainCoupler` stores the primary tissue domain and a secondary domain
and provides three hook points:

- `prepareSecondaryCoupling()`
  - runs before the upstream auxiliary domain advances
- `preparePrimaryCoupling()`
  - runs after the auxiliary-domain advance and before the primary myocardium
    solve
- `preparePostPrimaryCoupling()`
  - runs before a downstream domain such as ECG advances

This split supports staged domain updates without hard-wiring a single coupling
order into each domain.

## Other shipped coupling code

`heartBathInterfaceCoupler` is still compiled in this tree. It is bath-related
coupling code, but bath is not currently part of the active `core`
orchestration path.

## PVJ coupling family

`PVJMapper` owns the geometry work needed for 1D-to-3D exchange:

- locate myocardium cells associated with terminal Purkinje nodes
- gather tissue `Vm` at PVJ locations
- convert terminal currents into volumetric source terms
- distribute those source terms into the myocardium `sourceField`
- deposit terminal activation times into the myocardium activation field

`PVJCoupler` owns the family-level scaffolding:

- `networkCouplingEndpoint` discovery
- `couplingMode` parsing
- `PVJMapper`
- shared terminal coupling buffers

`ReactionDiffusionPvjCoupler` then applies the current resistive model at each PVJ:

```text
network terminal Vm  ----\
                          >  I_pvj = (V_network - V_tissue) / R_pvj
tissue Vm at PVJ     ----/
```

The coupler reuses internal buffers for:

- tissue PVJ voltages
- network PVJ voltages
- terminal coupling currents
- myocardium-side volumetric source equivalents

## Coupling modes

`ReactionDiffusionPvjCoupler` supports:

- `unidirectional`
  - one-way driving from the Purkinje network into the myocardium
  - the myocardium receives injected current
  - the network-side coupling buffers are zeroed before the network advance
- `bidirectional`
  - two-way exchange between network and myocardium
  - allows retrograde influence from tissue state through the resistive PVJ
    term

`EikonalPvjCoupler` supports:

- `unidirectional`
  - one-way transfer of Purkinje terminal activation times into the myocardium
    eikonal domain
- `bidirectional`
  - reserved, but currently stops with an "in development" message

## Runtime sequence

In the default staggered workflow:

1. The coupler samples the previous tissue and network states.
2. The conduction-system domain advances.
3. The coupler recomputes the PVJ exchange using the updated network state.
4. The myocardium receives the projected volumetric source and advances.

This is efficient and works well for one-way coupling. Bidirectional coupling
is stiffer and generally benefits from the iterative
`pimpleStaggeredElectrophysicsAdvanceScheme`; see
[../core/ARCHITECTURE.md](../core/ARCHITECTURE.md)
for the timestep-level tradeoffs.
