# Electrophysiology solvers

This directory contains OpenFOAM-based solvers for cardiac electrophysiology,
built on top of the `ionicModels` library.

```
solvers/
├── electroActivationFoam
├── eikonalElectroActivationFoam
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

This solver is intended for **1D/2D/3D tissue simulations** where full
electrophysiological dynamics are required.

---

## eikonalElectroActivationFoam

`eikonalElectroActivationFoam` is a **reduced-order electrophysiology solver**
for computing **activation times only**, based on an anisotropic
eikonal–diffusion formulation.

Instead of solving for the transmembrane voltage, the solver computes an
activation time field by solving a steady nonlinear eikonal–diffusion equation,
providing a computationally inexpensive approximation of the electrical wavefront
propagation.

Key features include:

- Anisotropic propagation via a conductivity (metric) tensor,
- Support for fibre-based conduction,
- Prescribed activation times on selected regions or patches,
- An optional stabilised Picard (defect-correction) formulation to improve
  convergence of the nonlinear solve.

This solver is intended for:

- Rapid estimation of activation times,
- Large-scale or parameter-sweep studies,
- Providing activation-time input to downstream electrophysiology or
  electrocardiographic models.

It does **not** model transmembrane voltage dynamics or ionic currents.

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

The three solvers address **different levels of electrophysiological modelling**:

- `electroActivationFoam`  
  Full tissue-scale electrophysiology (monodomain PDE + ionic ODEs).

- `eikonalElectroActivationFoam`  
  Reduced-order model computing activation times only, without voltage or ionic
  dynamics.

- `singleCellElectroActivationFoam`  
  Zero-dimensional single-cell simulations for ionic-model development and
  testing.

`singleCellElectroActivationFoam` can be viewed as a special case of
`electroActivationFoam` in which spatial coupling is removed. In the future,
`electroActivationFoam` may support a native single-cell mode, in which case
`singleCellElectroActivationFoam` may be deprecated.

---

## Notes

- `electroActivationFoam` and `singleCellElectroActivationFoam` rely on the
  `ionicModels` library for all ionic-model logic.
- `eikonalElectroActivationFoam` does not use ionic models, but shares geometry,
  material parameter, and stimulus concepts with the full solver.
- Selection of models, parameters, and numerical options is done via dictionaries
  at run-time.
- All solvers follow standard OpenFOAM build and execution patterns.
