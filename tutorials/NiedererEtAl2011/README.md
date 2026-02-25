# electroActivationFoam tutorial: Niederer et al. (2011) slab benchmark

This tutorial demonstrates the use of **electroActivationFoam** to reproduce the
cardiac tissue electrophysiology benchmark proposed by **Niederer et al. (2011)**,
which is widely used to verify cardiac electrophysiology simulators.

The benchmark is described in:

> S. A. Niederer et al.,  
> *Verification of cardiac tissue electrophysiology simulators using an N-version benchmark*,  
> Philosophical Transactions of the Royal Society A, 369:4331–4351, 2011.

---

## Benchmark overview

The Niederer et al. benchmark defines a **monodomain electrophysiology problem**
on a simple slab geometry with well-specified material parameters, stimulation
protocol, and evaluation metrics. The purpose is **code verification**, not model
validation.

Key characteristics of the benchmark:

- **Geometry:** rectangular slab (20 × 7 × 3 mm)
- **Model:** monodomain reaction–diffusion equation
- **Conductivity:** transversely isotropic
- **Cell model:** human ventricular epicardial cell
- **Stimulus:** localised cuboid stimulus applied at one corner
- **Metric:** activation time (first crossing of 0 mV)

---

## Case structure

```
NiedererEtAl2011/
├── 0
│   └── Vm
├── constant
│   ├── cardiacProperties
│   ├── timeIntegrationProperties
│   └── stimulusProtocol
├── system
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── setupNiedererEtAl2012
│   └── post-processing and reference scripts
├── Allrun
└── Allclean
```

---

## Solver options

This case can be run using either:

- **`electroActivationFoam`**  
  A full **monodomain reaction–diffusion solver**, reproducing the Niederer et al.
  benchmark as originally defined. This option is **more accurate** but
  computationally more expensive.

- **`eikonalElectroActivationFoam`**  
  A simplified **eikonal-based solver** that computes activation times directly
  by solving an anisotropic eikonal–diffusion equation. This option is
  **significantly cheaper**, but neglects detailed transmembrane dynamics and is
  intended for rapid estimation of activation times rather than full
  electrophysiological modelling.

The eikonal solver uses the same geometry, material parameters, and stimulus
location, but replaces the monodomain PDE with a steady nonlinear eikonal
formulation.

---

## Running the tutorial

By default, the tutorial runs the **monodomain solver** in serial:

```bash
./Allrun
```

To run the **eikonal solver** instead:

```bash
./Allrun eikonal
```

To run either solver in parallel:

```bash
./Allrun parallel
```

or, for the eikonal solver in parallel:

```bash
./Allrun eikonal parallel
```

---

## Post-processing and verification

Results can be post-processed to extract activation times and compared against
published benchmark data or other simulation codes.

---

## Purpose

This tutorial is intended for solver verification, regression testing, and
educational use, rather than physiological prediction.
