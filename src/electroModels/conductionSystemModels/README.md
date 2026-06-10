# conductionSystemModels

This directory contains runtime-selectable solver kernels used by
`conductionSystemDomain`. The domain owns the graph topology, PVJ metadata,
ionic-model state, and terminal-coupling buffers; the classes here implement
the numerical update applied to that graph state.

## Current contents

```text
src/electroModels/conductionSystemModels/
├── monodomain1DSolver/   # 1D cable-equation graph solver
├── eikonalSolver1D/      # Activation-time graph solver
├── restitutionEikonalSolver1D/  # Restitution-aware activation solver
└── README.md
```

## Available solvers

- `monodomain1DSolver`
  - Registered as `monodomain1DSolver`.
  - Advances graph-backed Purkinje state using ionic reaction terms plus an
    implicit cable-equation diffusion step.
  - The diffusion solve assumes a tree topology prepared by `conductionGraph`.

- `eikonalSolver1D`
  - Registered as `eikonalSolver1D`.
  - Computes nodal activation times only, using edge lengths and a prescribed
    wave speed `c0`.
  - Supports reduced-order conduction studies without full ionic state.

- `restitutionEikonalSolver1D`
  - Registered as `restitutionEikonalSolver1D`.
  - Tracks refractory state, recovery time, diastolic interval, action
    potential duration, and per-node activation history across timesteps.
  - Uses tabulated APD(DI) and CV(DI) restitution curves through
    `restitutionModel`.
  - Reports diagnostic fields for block, wavebreak, short-DI events, and
    minimum DI.

## Relationship to graph/domain code

The solver kernels depend on state owned by `conductionSystemDomain`, including:

- `conductionGraph`
- nodal `Vm1D`, `Iion1D`, `activationTime`, and solver diagnostic fields
- PVJ coupling currents prepared by electro couplers
- the runtime-selected ionic model for the graph nodes

That separation keeps graph construction and coupling bookkeeping in the domain
while leaving timestep logic in the solver implementations. See
[../electroDomains/README.md](../electroDomains/README.md) for the domain-level
API and [../electroCouplers/README.md](../electroCouplers/README.md) for PVJ
coupling.
