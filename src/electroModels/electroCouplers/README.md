# electroCouplers

The couplers that pass data between electrical domains, and the interfaces the domains expose to them.

## What's available

- Coupling endpoints (`electroDomainCouplingEndpoints.H`): `tissueCouplingEndpoint`, implemented by the myocardium; `networkCouplingEndpoint`, implemented by the Purkinje network; and `bathCouplingEndpoint`.
- `electroDomainCoupler`: the base class for couplers.
- The Purkinje–muscle junction (PVJ) family, which links the Purkinje network to the myocardium:
  - `pvjCoupler`, the family's base class, and `pvjMapper`, which finds the tissue cells at each junction and moves sources between the network and the tissue mesh;
  - `eikonalPvjCoupler`, `eikonalMonodomainPvjCoupler` and `reactionDiffusionPvjCoupler`. The last passes a resistive current,
    `(V_network − V_tissue) / R_pvj` in amperes with `R_pvj` in ohms, at each junction. In `couplingMode bidirectional` the
    network node loses that current from its volume `π ρ² L`, where `ρ` is `purkinjeFibreRadius` (default 20 µm) and `L` the node's control
    length, and the tissue receives the same current, so charge is conserved by construction. The
    [architecture notes](../ARCHITECTURE.md) give the equations.
  - `pvjDelay` (s, default 0): the anterograde junction delay. `eikonalPvjCoupler` seeds the tissue at a junction at the
    terminal's activation time plus `pvjDelay`; `eikonalMonodomainPvjCoupler` starts the Purkinje action-potential
    template that much later. `reactionDiffusionPvjCoupler` refuses a non-zero value: its delay comes from `rPvj` and
    the junction load. Measured delays are about 3–13 ms (rabbit, dog).
  - `seedProfile` (`eikonalPvjCoupler`, default `front`): the seed time of the support
    cells. `flat` holds them at the terminal time; `front` adds `rho/c0` and `curvedFront`
    the source-free radial solution, with `rho` the metric distance to the terminal. A
    time away from the terminal depends on `pvjRadius` linearly for `flat`, through a
    logarithm for `front`, and not at all for `curvedFront`.

The bath has no coupler: `extracellularPotentialDomain` shares `phiE` with the myocardium directly.

## Folders

```text
src/electroModels/electroCouplers/
├── electroDomainCouplingEndpoints.H
├── electroDomainCoupler.{H,C}
└── pvjCoupler/
    ├── pvjCoupler.{H,C}
    ├── pvjMapper.{H,C}
    ├── eikonal/
    ├── eikonalMonodomain/
    └── reactionDiffusion/
```
