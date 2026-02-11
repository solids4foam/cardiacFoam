# cardiacFoam

**cardiacFoam** is an open-source OpenFOAM toolbox for cardiac electro-mechanical
simulation, developed as part of the **XenoSim** project.

It is designed as a cardiac-focused extension of
[solids4foam](https://github.com/solids4foam/solids4foam), building directly on
its finite-volume solid mechanics and fluid–structure interaction (FSI)
capabilities, and extending them with cardiac electrophysiology and
electromechanical coupling models.

The code is under active development and should be considered **work in
progress**.

---

## Scope and objectives

The long-term aim of **cardiacFoam** is to provide a unified finite-volume
framework for coupled cardiac electrophysiology, solid mechanics, and
fluid–structure interaction within OpenFOAM.

Specifically, the toolbox targets:

- Cardiac electrophysiology modelling (reaction–diffusion and reduced models)
- Electromechanical coupling for active myocardium
- Integration with FSI formulations available in solids4foam
- Research-grade, extensible solvers suitable for whole-heart simulations

---

## What does it currently contain?

The toolbox currently includes:

- **Electrophysiology solvers**
  - Monodomain reaction–diffusion model
  - Eikonal–diffusion model
  - Minimal Bueno-Orovio ionic model

- **Tutorial cases**
  - Niederer et al. (2011) benchmark slab test case
  - Additional example cases demonstrating electro-only workflows

- **Ongoing development**
  - Electromechanical coupling formulations
  - Cardiac-specific constitutive and activation models
  - Closer integration with solids4foam FSI capabilities

---

## Dependencies and compatibility

- **OpenFOAM:**
  Compatible with **OpenFOAM-v2312 through OpenFOAM-v2512**

- **solids4foam:**
  Required as a dependency and intended to be included as a **Git submodule**.
  cardiacFoam builds directly on the solid mechanics and FSI infrastructure
  provided by solids4foam.

- **PETSc:**
  Used for nonlinear and coupled solution strategies where required.

---

## Installation

1. Source a supported OpenFOAM version (v2312–v2512).
2. Clone the repository and initialise submodules:

   ```bash
   git clone https://github.com/<organisation>/cardiacFoam.git
   cd cardiacFoam
   git submodule update --init --recursive
   ```

3. Compile the toolbox:

   ```bash
   ./Allwmake
   ```

---

## Running a tutorial

To run the Niederer benchmark slab case:

```bash
cd tutorials/NiedererBenchmarkSlab
./Allrun
```

---

## Status

**cardiacFoam is research software under active development.**
APIs, solvers, and case setups may change as the project evolves.

Contributions, feedback, and collaboration are very welcome.

---

## Contact

For questions or collaboration inquiries, please email [philip.cardiff@ucd.ie](mailto:philip.cardiff@ucd.ie)
