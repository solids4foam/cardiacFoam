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

## Cardiac model and material parameters

The file `constant/cardiacProperties` defines the electrophysiology model and
material parameters consistent with Niederer et al. (2011), including anisotropic
conductivity, membrane capacitance, surface-to-volume ratio, ionic model choice,
and tissue type.

---

## Time integrations settings

The file `constant/timeIntegrationProperties` controls the numerical strategy,
including explicit/implicit coupling, CFL limits, ODE solver selection, and
integration tolerances.

---

## Stimulation protocol

The file `constant/stimulusProtocol` defines a localised cuboid stimulus applied
at one corner of the slab for 2 ms, matching the benchmark definition.

---

## Running the tutorial

```bash
./Allrun
```

To clean the case:

```bash
./Allclean
```

---

## Post-processing and verification

Results can be post-processed to extract activation times and compared against
published benchmark data or other simulation codes.

---

## Purpose

This tutorial is intended for solver verification, regression testing, and
educational use, rather than physiological prediction.
