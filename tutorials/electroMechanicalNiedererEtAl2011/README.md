# NiedererEtAl2012_EM: Electro-Mechanical Coupling Tutorial

## Purpose

This tutorial demonstrates the electro-mechanical coupling framework in
cardiacFoam using the `electroMechanicalModel` physics model. It extends the
electro-only NiedererEtAl2012 benchmark to include a solid mechanics region
coupled to the electrophysiology region.

The case uses a 20x3x7 mm tissue slab (Niederer et al. 2011 benchmark
geometry) with:

- **Electro region**: monodomain reaction-diffusion PDE with the TNNP ionic
  model
- **Solid region**: nonlinear total Lagrangian solid solver with the
  `electroMechanicalLaw` (neo-Hookean passive + active tension)

The two regions are coupled sequentially: after each electro solve, the
intracellular calcium concentration (Cai) is extracted from the ionic model
and converted to an active tension field (Ta) using a simple linear model.
This Ta field is passed to the solid region where the `electroMechanicalLaw`
adds it as a fibre-aligned active stress component.

## Running

```bash
./Allrun
```

## Case Structure

```text
constant/
    physicsProperties              # type electroMechanicalModel
    electroMechanicalProperties    # coupling scheme selection
    electro/                       # electro region dictionaries
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

## Coupling parameters

Set in `constant/electroMechanicalProperties`:

```
kTa             1e7;       # Pa per (Cai unit); TNNP Cai is in mM
CaiThreshold    0.0002;    # mM; TNNP resting Cai is ~0.0002 mM
```

Active tension: `Ta = kTa * max(Cai - CaiThreshold, 0)`.
This simple linear model is a placeholder for a dedicated active tension model.
