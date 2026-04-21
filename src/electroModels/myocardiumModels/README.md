# myocardiumModels

This directory contains myocardium-side numerical kernels and related electro
models. The domain-side ownership lives in `electroDomains/myocardiumDomain/`,
while the classes here provide the diffusion, reduced-order, or single-cell
solve implementations used by those workflows.

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
  - In the current architecture, the eikonal myocardium path is surfaced
    through `EikonalMyocardiumDomain`, not by forcing the classic
    reaction-diffusion `MyocardiumDomain` to own it directly.

## ODE-only workflow

- `SingleCellSolver`
  - Registered as `singleCellSolver`.
  - Advances one integration point with a runtime-selected ionic model and no
    spatial PDE.
  - Used for ionic-model testing, calibration, and waveform generation rather
    than tissue-scale propagation.

## Relationship to myocardium-domain code

`MyocardiumDomain` owns the reaction-diffusion tissue path, while
`EikonalMyocardiumDomain` owns the reduced-order activation-time path.

Together, the myocardium-domain layer owns:

- the tissue fields (`Vm`, `Iion`, `sourceField`, `activationTime`)
- the ionic model
- stimulus handling
- export/post-processing hooks

The solver classes in this directory own the numerical kernel associated with
the chosen tissue formulation or workflow. See
[../electroDomains/README.md](../electroDomains/README.md) for the domain-level
state and lifecycle.
