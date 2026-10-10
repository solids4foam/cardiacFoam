# numericsModels

The discretisation pieces the tissue solvers share: how a sealed heart wall and a heart-bath interface are treated on the mesh. They depend only on OpenFOAM's finite-volume library, so they carry no knowledge of monodomain, bidomain or bath; the solvers call them.

## What's available

| File | What it does |
|---|---|
| `sealedWall/insulatedFaceConductivity` | face conductivity with zero values on sealed walls, so no current crosses them |
| `sealedWall/conormalZeroFluxFvPatchScalarField` | boundary condition whose gradient sets the conormal flux of a potential to zero |
| `interface/extracellularFaceConductivity.H` | face conductivity across the heart-bath interface: linear, harmonic, distance-weighted harmonic and conormal harmonic (header-only) |
| `interface/exposedPhiETrace` | value of the extracellular potential on a heart-bath face from two-sided flux continuity, and its use as the heart's boundary trace |
| `interface/interfaceGradient` | cell gradient of a potential whose normal derivative jumps across the heart-bath interface |

## Folders

```text
src/numericsModels/
├── sealedWall/
└── interface/
```

## What this does not own

- The equations and their coupling: the solvers and domains in [electroModels](../electroModels/README.md).
