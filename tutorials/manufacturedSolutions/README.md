# Manufactured Solutions Tutorials

This folder groups the manufactured-solution verification cases by model scope.

## Cases

- `monodomainPseudoECG` : spatial manufactured-solution verification for the monodomain solver with pseudo-ECG output
- `bidomain` : spatial manufactured-solution verification for the bidomain solver
- `singleCellMonodomain` : single-cell manufactured-solution verification for the monodomain ionic workflow
- `singleCellBidomain` : single-cell manufactured-solution verification for the bidomain ionic workflow

## Naming Pattern

The subdirectory names follow:

- `<model>` for spatial PDE manufactured cases
- `singleCell<model>` for zero-dimensional manufactured cases

This keeps all manufactured cases under one parent while preserving the distinction between PDE and single-cell workflows.

## Driver Usage

Run the whole manufactured suite:

```bash
./Allrun
```

Run only selected cases:

```bash
./Allrun monodomainPseudoECG bidomain
./Allrun singleCellMonodomain singleCellBidomain
```
