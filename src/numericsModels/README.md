# numericsModels

The discretisation pieces the tissue solvers share: how a sealed heart wall and a heart-bath interface are treated on the mesh. They depend only on OpenFOAM's finite-volume library, so they carry no knowledge of monodomain, bidomain or bath; the solvers call them.

## What's available

| File | What it does |
|---|---|
| `sealedWall/sealedWallPolicy` | which patches of a potential are sealed, and the per-patch wall-gradient states (stock, sealed, value) a field's owner publishes |
| `sealedWall/insulatedFaceConductivity` | face conductivity with zero values on sealed walls, so no current crosses them |
| `sealedWall/sealedWallLeastSquaresGrad` | gradient scheme `sealedWallLeastSquares`: least squares on the point-neighbour stencil without sealed-face values in wall cells, `leastSquares` elsewhere |
| `sealedWall/sealedWallLeastSquaresVectors` | the wall-cell least-squares vectors of one field, cached on its mesh |
| `interface/extracellularFaceConductivity.H` | face conductivity across the heart-bath interface: linear, harmonic, distance-weighted harmonic and conormal harmonic (header-only) |
| `interface/exposedPhiETrace` | value of the extracellular potential on a heart-bath face from two-sided flux continuity, and its use as the heart's boundary trace |
| `interface/interfaceGradient` | cell gradient of a potential whose normal derivative jumps across the heart-bath interface |
| `pureNeumann/compatibilityImbalance` | reports the discrete imbalance of a pure-Neumann potential equation, which its reference cell absorbs; the equation is not changed |
| `probedOperator/probedOperator` | a linear cell-to-cell operator stored as a sparse matrix by graph-coloured probing of the function that applies it, with its transpose |
| `diffusion/conormalStabilisation` | term added to a tensor Laplacian that makes the face-conormal magnitude its implicit coefficient, so the correctors converge when the conductivity is oblique to the faces; zero at convergence |
| `advection/facePeclet` | face Peclet number of an advection-diffusion operator with a tensor face diffusivity, and its log line |

## Folders

```text
src/numericsModels/
├── sealedWall/
├── interface/
├── pureNeumann/
├── diffusion/
├── probedOperator/
└── advection/
```

## What this does not own

- The equations and their coupling: the solvers and domains in [electroModels](../electroModels/README.md).
