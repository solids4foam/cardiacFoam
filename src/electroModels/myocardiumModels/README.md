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
├── singleCellSolver/   # ODE-only single-cell workflow
└── README.md
```

## Tissue-domain solvers

- `monodomainSolver`
  - Registered as `monodomainSolver`.
  - Solves the standard tissue reaction-diffusion problem for `Vm`.
  - Owns the monodomain conductivity tensor used by explicit and implicit
    diffusion updates.

- `bidomainSolver`
  - Registered as `bidomainSolver`.
  - Solves the coupled tissue system for `Vm` and extracellular potential
    `phiE`.
  - Owns intracellular and extracellular conductivity tensors and the `phiE`
    field exposed through `myocardiumDomain`.

There is no eikonal kernel in this directory. The dictionary entry
`myocardiumSolver eikonalSolver` selects `electrophysiologyModel`, which builds
`eikonalMyocardiumDomain` (in `../electroDomains/myocardiumDomain/`). That domain
assembles the eikonal-diffusion equation itself; there is no separate eikonal
`myocardiumSolver`.

## ODE-only workflow

- `singleCellSolver`
  - Registered as `singleCellSolver` in the parent `electroModel` table, not in
    the `myocardiumSolver` diffusion-solver table.
  - Advances one integration point with a runtime-selected ionic model and no
    spatial PDE.
  - Used for ionic-model testing, calibration, and waveform generation rather
    than tissue-scale propagation.

The solver classes in this directory own the numerical kernel for the chosen
tissue formulation. See
[../electroDomains/README.md](../electroDomains/README.md) for the domain-level
state and lifecycle.
