# cardiacFoam tutorial: Niederer et al. (2011) slab benchmark

This tutorial demonstrates the use of **cardiacFoam** to reproduce the
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
NiedererEtAl2012/
├── 0
│   └── Vm
├── constant
│   ├── electroProperties
│   ├── physicsProperties
│   ├── activeTensionProperties
│   └── polyMesh
├── system
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── setupNiedererEtAl2012
│   └── automated simulation pipeline
│   └── post-processing and reference scripts
├── Allrun
└── Allclean
```

---

## Solver options

This case runs with **`cardiacFoam`**.

---

## Running the tutorial

By default, the tutorial runs in serial:

```bash
./Allrun
```

To run in parallel:

```bash
./Allrun parallel
```

or, for the eikonal solver in parallel:

---

## Post-processing and verification

Results can be post-processed to extract activation times and compared against
published benchmark data or other simulation codes.
The code presented in the setup runs several simulations and run the post-process
for mesh and time dependant activation times. 


---

## Purpose

This tutorial is intended for solver verification, regression testing, and
educational use, rather than physiological prediction.
