# modules/correctedFvPatchFields architecture

This module provides `fixedValueCorrected` and `fixedGradientCorrected` boundary
conditions (scalar and vector) so `cardiacFoam` can use them regardless of
whether it's linked against the full `solids4foam` or the lightweight
`modules/physicsModel` stand-in.

## Why this exists

`pointCellsLeastSquares` extends the least-squares gradient stencil through mesh
points (cell-point-cell connectivity), not just face-neighbors. Standard
OpenFOAM boundary conditions compute `snGrad()` with a plain orthogonal
approximation, which is inconsistent with that extended stencil and silently
drops accuracy at boundary cells (see solids4foam's own `patchTest` tutorial).
`solids4foam` already fixes this with a small, genuinely standalone correction
(`patchCorrectionVectors`) applied by a family of "corrected" boundary
conditions -- but that code lives inside the much larger `solids4FoamModels`
library (mechanical models, FSI, RBF mesh motion, ...), which this project only
links when `USE_LIGHTWEIGHT_PHYSICSMODEL=0`.

This module builds *just* the correction closure as its own tiny library, so it
is available unconditionally -- no dependency on which `physicsModel` backend
is in use.

## Directory structure

```text
modules/correctedFvPatchFields/
├── Make/files      -- sources itself directly from modules/solids4foam
│                      (not copied, so upstream fixes flow through on a
│                      normal 'git submodule update')
├── Make/options     -- only needs -lfiniteVolume
└── README.md
```

## What it provides

Runtime-selectable boundary condition types (usable in any `0/<field>` dict via
their `type` entry, no C++ changes needed elsewhere):

- `fixedValueCorrected` (scalar and vector) -- Dirichlet with non-orthogonal
  correction.
- `fixedGradientCorrected` (scalar and vector) -- Neumann with non-orthogonal
  correction.

## Scope

Deliberately narrow: only the boundary-condition correction closure, not any
other part of solids4foam's numerics. If more corrected-BC types are needed
later, add their source paths to `Make/files` the same way.
