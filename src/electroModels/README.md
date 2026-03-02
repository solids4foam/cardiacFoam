# electroModels

This directory contains the **electroModels library** for cardiac electrophysiology
in cardiacFoam. It replaces the earlier solver-based structure with a **modular,
run-time selectable model framework**, following the same design philosophy as
solids4foam (e.g. `solidModels`, `fluidModels`).

---

## Overview

Instead of separate solvers for different electrophysiology formulations,
cardiacFoam now provides a **single solver framework** with **run-time selectable
electroModels**.

Each electroModel encapsulates a specific formulation, such as:

- Monodomain reaction–diffusion
- Eikonal–diffusion activation
- Single-cell ionic dynamics

This approach improves:

- Code reuse
- Extensibility
- Consistency with solids4foam design patterns

---

## Library structure

```bash
electroModels/
├── electroModel/                 # Base class
├── monoDomainElectro/            # Monodomain PDE + ionic ODEs
├── eikonalDiffusionElectro/      # Eikonal–diffusion model
├── singleCellElectro/            # 0D ionic model driver
├── lnInclude/                   # OpenFOAM linking headers
└── Make/
```

---

## electroModel (base class)

`electroModel` defines the **common interface** for all electrophysiology models.

Responsibilities include:

- Managing fields (e.g. transmembrane voltage `Vm`)
- Handling time integration
- Interfacing with ionic models (where applicable)
- Providing a consistent API to the solver

All specific models derive from this base class and are selected at run-time.

---

## monoDomainElectro

Implements the **monodomain reaction–diffusion equation**:

- Solves for transmembrane voltage `Vm`
- Couples to ionic models (`ionicModels` library)
- Supports tissue-scale simulations (1D/2D/3D)

This corresponds to the functionality previously provided by
`electroActivationFoam`.

---

## eikonalDiffusionElectro

Implements a **reduced-order eikonal–diffusion model**:

- Computes activation times only
- Does not solve for voltage or ionic currents
- Suitable for fast propagation simulations and large-scale studies

---

## singleCellElectro

Implements a **zero-dimensional (0D) single-cell model**:

- Solves ionic model ODEs only
- No spatial coupling (no PDE)
- Used for testing and validating ionic models

This corresponds to the previous `singleCellElectroActivationFoam`.

---

## Run-time selection

Electrophysiology models are selected at run-time via dictionaries, e.g.:

```plaintext
electroModel   monoDomainElectro;
```

This mirrors the solids4foam approach:

```plaintext
solidModel     linearElastic;
fluidModel     incompressible;
```

---

## Design philosophy

The electroModels framework:

- Follows the **run-time selection mechanism** used in OpenFOAM
- Mirrors the **modular physics design of solids4foam**
- Enables easy addition of new electrophysiology formulations
- Avoids duplication of solver infrastructure

---

## Relationship to ionicModels

- `monoDomainElectro` and `singleCellElectro` use the `ionicModels` library
- `eikonalDiffusionElectro` does not require ionic models

---

## Notes

- The solver layer is now thin and delegates physics to electroModels
- Future extensions (e.g. bidomain) can be implemented as additional
  electroModels
- This structure enables tighter coupling with solids4foam in future

---

## Summary

| Model                  | Description                                 |
|------------------------|---------------------------------------------|
| monoDomainElectro      | Full PDE + ionic ODEs                       |
| eikonalDiffusionElectro| Activation times only                       |
| singleCellElectro      | 0D ionic model                              |

---

This refactoring aligns cardiacFoam with modern OpenFOAM and solids4foam
design practices, improving maintainability and extensibility.
