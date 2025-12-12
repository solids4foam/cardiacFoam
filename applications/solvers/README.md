# Electrophysiology solvers

This directory contains OpenFOAM-based solvers for cardiac electrophysiology,
built on top of the `ionicModels` library.

```
solvers/
├── electroActivationFoam
└── singleCellElectroActivationFoam
```

---

## electroActivationFoam

`electroActivationFoam` is the main electrophysiology solver.
It solves a **coupled PDE–ODE system**:

- **PDE:** reaction–diffusion equation for the transmembrane voltage `Vm`
  (e.g. monodomain formulation),
- **ODEs:** ionic cell model equations, provided at run-time via
  `Foam::ionicModel`.

Key features include:

- Run-time selectable ionic models,
- Explicit and implicit time-integration paths,
- Support for manufactured solutions (verification),
- Activation-time tracking.

This solver is intended for **1D/2D/3D tissue simulations**.

---

## singleCellElectroActivationFoam

`singleCellElectroActivationFoam` is a **single-cell driver** for ionic models.

In this case, the voltage equation is no longer a PDE; instead:

- The transmembrane voltage `Vm` is solved **inside the ionic model ODE system**,
- The solver advances a **single integration point** in time,
- Results are written as a time series using `ionicModelIO`.

This solver is primarily intended for:

- Testing and validating ionic models,
- Reproducing single-cell action potential traces,
- Debugging and development.

---

## Relationship between the solvers

`singleCellElectroActivationFoam` is effectively a **special case** of
`electroActivationFoam`:

- Both use the same `ionicModel` infrastructure,
- The difference is whether `Vm` is solved via a PDE or as part of the ODE system.

In the future, `electroActivationFoam` is expected to support a single-cell mode
directly. When this is implemented, `singleCellElectroActivationFoam` will likely
be deprecated and removed.

---

## Notes

- Both solvers rely on the `ionicModels` library for all ionic-model logic.
- Selection of the ionic model and stimulus protocol is done via dictionaries
  at run-time.
- The solvers follow standard OpenFOAM build and execution patterns.
