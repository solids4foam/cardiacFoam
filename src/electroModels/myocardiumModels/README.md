# myocardiumModels

The solvers for the heart tissue, plus the single-cell solver.

## What's available

- `monodomainSolver`: reaction–diffusion for `Vm`.
- `bidomainSolver`: the coupled system for `Vm` and the extracellular potential `phiE`.
- `singleCellSolver`: one cell and no tissue, running an ionic model on its own. It can also compute active tension.

The eikonal tissue model has no solver here: `eikonalMyocardiumDomain`, in [electroDomains](../electroDomains/README.md), solves it itself.

## Folders

```text
src/electroModels/myocardiumModels/
├── monodomainSolver/
├── bidomainSolver/
└── singleCellSolver/
```
