# cardiacFoam

**cardiacFoam** is an open-source OpenFOAM toolbox for cardiac electro-mechanical
simulation, developed as part of the **XenoSim** project.

It is designed as a cardiac-focused extension of
[solids4foam](https://github.com/solids4foam/solids4foam), building on its
finite-volume solid mechanics and fluid–structure interaction (FSI)
capabilities, while extending them with cardiac electrophysiology and
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

## Design philosophy

**cardiacFoam** follows the same modular, run-time selectable design philosophy
as solids4foam.

- Electrophysiology models are implemented as run-time selectable
  **`electroModels`**
- This mirrors the **`solidModels`** and **`fluidModels`** architecture in
  solids4foam
- New electrophysiology formulations (e.g. monodomain, eikonal, etc.) can be
  added in a consistent and extensible way

This approach enables flexible model selection at run-time and promotes
code reuse and extensibility across different cardiac modelling workflows.

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

## Dependencies and usage modes

- **OpenFOAM:**
  Compatible with **OpenFOAM-v2312 through OpenFOAM-v2512**

- **PETSc (optional):**
  Used for nonlinear and coupled solution strategies in solids4foam.

### solids4foam: optional dependency

**cardiacFoam can be used in two modes:**

#### 1. With full solids4foam (required for electromechanics / FSI)

If a full solids4foam installation is available, cardiacFoam will use it to
access:

- Solid mechanics models
- Fluid–structure interaction (FSI)
- Advanced constitutive laws and coupling infrastructure

This is required for:

- Electromechanical simulations
- Fully coupled heart mechanics + flow problems

You can specify the installation via:

```bash
export SOLIDS4FOAM_INST_DIR=/path/to/solids4foam
```

If not explicitly set, cardiacFoam will attempt to detect a solids4foam
installation (e.g. via submodules).

---

#### 2. Without solids4foam (lightweight electrophysiology-only mode)

If solids4foam is not available, cardiacFoam automatically falls back to a
lightweight internal **`physicsModel`** implementation.

This mode supports:

- Standalone electrophysiology simulations
- Monodomain and eikonal models
- Ionic models and activation dynamics

but **does not include**:

- Solid mechanics
- Electromechanical coupling
- FSI capabilities

---

### Automatic dependency resolution

Dependency handling is managed by:

```bash
etc/resolveSolids4Foam.sh
```

This script:

- Detects whether a valid `SOLIDS4FOAM_INST_DIR` is available
- Falls back to the internal `modules/physicsModel` if not
- Sets up the required include paths and libraries accordingly

This allows the same codebase to seamlessly support both:

- full multi-physics workflows (with solids4foam), and
- lightweight electro-only simulations (without it)

---

## Installation

1. Source a supported OpenFOAM version (v2312–v2512)

2. Clone the repository:

    ```bash
    git clone https://github.com/<organisation>/cardiacFoam.git
    cd cardiacFoam
    ```

3. (Optional) Initialise submodules (recommended for full functionality):

    ```bash
    git submodule update --init --recursive
    ```

4. Compile the toolbox:

    ```bash
    ./Allwmake
    ```

    During compilation, the appropriate solids4foam or lightweight backend will be
    automatically selected.

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

For questions or collaboration inquiries, please email
[philip.cardiff@ucd.ie]
