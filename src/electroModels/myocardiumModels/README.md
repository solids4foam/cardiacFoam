# myocardiumModels

This directory contains runtime-selectable numerical kernels used by
`MyocardiumDomain`. The domain owns the tissue fields and lifecycle; the
classes here provide the diffusion or reduced-order solve applied to those
fields.

## Current contents

```text
src/electroModels/myocardiumModels/
├── monodomainSolver/   # Single-potential tissue PDE
├── bidomainSolver/     # Coupled Vm and phiE tissue PDE
├── eikonalSolver/      # Reduced-order activation-time model
├── singleCellSolver/   # ODE-only single-cell workflow
└── README.md
```

## Tissue-domain solvers

- `MonodomainSolver`
  - Registered as `monodomainSolver`.
  - Solves the standard tissue reaction-diffusion problem for `Vm`.
  - Owns the monodomain conductivity tensor used by explicit and implicit
    diffusion updates.

- `BidomainSolver`
  - Registered as `bidomainSolver`.
  - Solves the coupled tissue system for `Vm` and extracellular potential
    `phiE`.
  - Owns intracellular and extracellular conductivity tensors and the `phiE`
    field exposed through `MyocardiumDomain`.

- `EikonalSolver`
  - Registered as `eikonalSolver`.
  - Computes activation times with a reduced-order anisotropic eikonal
    formulation rather than a full ionic-PDE solve.
  - Useful for fast propagation studies where only activation timing is needed.

## ODE-only workflow

- `SingleCellSolver`
  - Registered as `singleCellSolver`.
  - Advances one integration point with a runtime-selected ionic model and no
    spatial PDE.
  - Used for ionic-model testing, calibration, and waveform generation rather
    than tissue-scale propagation.

## Relationship to `MyocardiumDomain`

`MyocardiumDomain` owns:

- the tissue fields (`Vm`, `Iion`, `sourceField`, `activationTime`)
- the ionic model
- stimulus handling
- export/post-processing hooks

The solver classes in this directory own only the numerical kernel associated
with the chosen tissue formulation. See
[../electroDomains/README.md](../electroDomains/README.md) for the domain-level
state and lifecycle.
