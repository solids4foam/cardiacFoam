# Fiber-Aligned Anisotropic Mesh Pipeline

This tutorial demonstrates the generation and simulation of fiber-aligned anisotropic polyhedral dual meshes for cardiac electrophysiology using the Finite Volume Method (FVM) solver `cardiacFoam`.

## Directory Structure
- **benchmarks/**: Validation and reference cases (e.g., Niederer et al. 2011/2012 benchmark).
- **meshes/**: Raw geometry files, Gmsh output (`.msh`), background metrics (`.pos`), and cobiveco input/output (`.vtu`).
- **scripts/**: Python scripts for mesh generation and bash scripts for automating the pipeline.
- **simulations/**: Active simulation cases for comparing tetrahedral and polyhedral meshes at different resolutions.
- **templates/**: Base template cases used to initialize the simulation folders.

## Mesh scaling and content generation

The UVC coordinates, fiber/sheet/normal directions, and the conductivity tensor
are **scale-invariant**: none of the generators reads mesh geometry (cell
volumes, face areas, `deltaCoeffs`, or absolute lengths) into the content.

- Transmural/apicobasal coordinates are clamped to `[0,1]` and built from
  distance *ratios* + normalized gradients (`cobivecco-OpenFOam/src/fiberFoam/coordinates/`).
- The LDRBM fiber model consumes only `tm in [0,1]` and normalized frame
  vectors, emitting unit vectors (`.../fibers/ldrbmFiberModel.C`).
- MATLAB Cobiveco emits the same normalized `[0,1]` UVC.
- `setCardiacConductivity` builds `df*(f⊗f) + ds*(s⊗s) + dn*(n⊗n)` from physical
  constants in `cardiacCoreDict` (`cardiacCoreStandalone/src/setCardiacConductivity/`).

`transformPoints` moves only `constant/polyMesh/points`; it never touches field
values or cell ordering. Therefore **generating content and then scaling the mesh
produces identical content to scaling first** — content never needs regenerating
after a scale.

**Canonical scale:** both pipelines scale the native `~4x5x4`-unit ellipsoid by
`(0.02 0.02 0.02)` to `~8x10x8 cm`, applied once as the final step after content
generation. (A previous `0.001` factor in the tet path produced a `~4x5x4 mm`
heart — 20x too small — and was a mid-test artifact, not an intended value.)

See `docs/superpowers/specs/2026-07-22-fiber-mesh-scale-invariance-design.md` for
the full source audit.
