# NiedererEtAl2011_EM: Electro-Mechanical Coupling Tutorial

## Purpose

This tutorial demonstrates the electro-mechanical coupling framework in
cardiacFoam using the `electroMechanicalModel` physics model. It extends the
electro-only NiedererEtAl2011 benchmark to include a solid mechanics region
coupled to the electrophysiology region.

Unlike the electro-only tutorials, this case requires a full `solids4foam`
build of `cardiacFoam`. It does not run in lightweight mode, because
`libelectroMechanicalModels` is only compiled when `cardiacFoam` is built
against a compiled `solids4foam` installation.

The case uses a 20x3x7 mm tissue slab (Niederer et al. 2011 benchmark
geometry) with:

- **Electro region**: monodomain reaction-diffusion PDE with the TNNP ionic
  model, selected through `myocardiumSolver monodomainSolver`
- **Solid region**: nonlinear total Lagrangian solid solver with the
  `electroMechanicalLaw` (neo-Hookean passive + active tension)

The two regions are coupled sequentially: after each electro solve, the
intracellular calcium concentration (Cai) is extracted from the ionic model
and converted to an active tension field (Ta) using a simple linear model.
This Ta field is passed to the solid region where the `electroMechanicalLaw`
adds it as a fibre-aligned active stress component.

## Running

Build requirement:

- compile `solids4foam`
- rebuild `cardiacFoam` in full mode (`./Allwmake` without the lightweight
  fallback)

Then run:

```bash
./Allrun
```

## Case Structure

```text
constant/
    physicsProperties              # type electroMechanicalModel
    electroMechanicalProperties    # coupling scheme selection
    electro/electroProperties      # current myocardiumSolver electro region dictionary
    solid/                         # solid region dictionaries
system/
    controlDict                    # shared time control
    blockMeshDict                  # shared mesh definition
    electro/                       # electro fvSchemes/fvSolution
    solid/                         # solid fvSchemes/fvSolution
0/
    solid/                         # solid initial conditions (D, f0, f0f)
```

The electro region creates its fields internally (Vm, etc.) with default
values. The mesh is generated once via `blockMesh` and copied to both regions.

## Electro dictionary

The electro region uses the same canonical dictionary keys as the electro-only
Niederer case:

```cpp
myocardiumSolver monodomainSolver;

monodomainSolverCoeffs
{
    ionicModel TNNP;
    tissue epicardialCells;
    solutionAlgorithm explicit;

    externalStimulus { ... }
}
```

## Coupling parameters

Set in `constant/electroMechanicalProperties`:

```text
kTa             1e7;       # Pa per (Cai unit); TNNP Cai is in mM
CaiThreshold    0.0002;    # mM; TNNP resting Cai is ~0.0002 mM
```

Active tension: `Ta = kTa * max(Cai - CaiThreshold, 0)`.
This simple linear model is a placeholder for a dedicated active tension model.
