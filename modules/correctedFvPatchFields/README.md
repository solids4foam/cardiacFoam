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
├── src/numerics/patchCorrectionVectors/    -- vendored copy from solids4foam
├── src/numerics/compatibilityFunctions/    -- vendored copy from solids4foam
├── src/fvPatchFields/fixedValueCorrected/  -- vendored copy from solids4foam
├── src/fvPatchFields/fixedGradientCorrected/ -- vendored copy from solids4foam
├── Make/files
├── Make/options     -- only needs -lfiniteVolume -lmeshTools
└── README.md
```

These five files are **vendored copies**, not sourced by relative path from
`modules/solids4foam`. That submodule is not checked out in CI (all three
workflows use `submodules: false`) and is also absent under
`USE_LIGHTWEIGHT_PHYSICSMODEL=1`, so a build-time reference to
`../solids4foam/...` fails in both cases -- this library needs to be buildable
unconditionally, which for these five self-contained files (no dependency on
the rest of `solids4FoamModels`) means owning a copy rather than reaching
across a submodule boundary that may not exist. If solids4foam's upstream
copies of these files change, re-copy them here manually; there is no
automatic sync.

`fixedValueCorrectedFvPatchScalarField.C`'s `write()` here also carries a bug
fix (calls `fixedValueFvPatchScalarField::write()` instead of
`fvPatchField<scalar>::write()`, matching how every other `fixedValue`-derived
BC serialises its `value` entry -- the original dropped it, breaking
`reconstructPar` after any parallel run using this BC) that has not been
applied to the `modules/solids4foam` submodule itself.

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
