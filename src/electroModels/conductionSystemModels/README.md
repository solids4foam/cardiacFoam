# conductionSystemModels

This directory contains runtime-selectable solver kernels used by
`ConductionSystemDomain`. The domain owns the graph topology, PVJ metadata,
ionic-model state, and terminal-coupling buffers; the classes here implement
the numerical update applied to that graph state.

## Current contents

```text
src/electroModels/conductionSystemModels/
├── monodomain1DSolver/   # 1D cable-equation graph solver
├── eikonalSolver1D/      # Activation-time graph solver
└── README.md
```

## Available solvers

- `Monodomain1DSolver`
  - Registered as `monodomain1DSolver`.
  - Advances graph-backed Purkinje state using ionic reaction terms plus an
    implicit cable-equation diffusion step.
  - In the current nested implementation, the diffusion solve assumes a tree
    topology prepared by `conductionGraph` and uses graph traversal data built
    by the domain.

- `EikonalSolver1D`
  - Registered as `eikonalSolver`.
  - Computes nodal activation times only, using edge lengths and a prescribed
    wave speed `c0`.
  - Intended for reduced-order conduction studies where full ionic state is
    unnecessary.

## Relationship to graph/domain code

The solver kernels depend on state owned by `ConductionSystemDomain`, including:

- `conductionGraph`
- nodal `Vm1D`, `Iion1D`, and `activationTime`
- PVJ coupling currents prepared by electro couplers
- the runtime-selected ionic model for the graph nodes

That separation keeps graph construction and coupling bookkeeping in the domain
while leaving timestep logic in the solver implementations. See
[../electroDomains/README.md](../electroDomains/README.md) for the domain-level
API and [../electroCouplers/README.md](../electroCouplers/README.md) for PVJ
coupling.
