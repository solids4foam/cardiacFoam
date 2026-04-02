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
  `electroMechanicalLaw` (Guccione passive + active tension)

Currently, the two regions are solved sequentially with no field exchange (Stage
1 shell coupling). Active tension coupling will be added in a future stage.

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
