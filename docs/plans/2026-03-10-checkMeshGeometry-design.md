# Design: checkMeshGeometry utility

Date: 2026-03-10

## Purpose

A standalone OpenFOAM utility that reads an existing `constant/polyMesh`, detects
whether the mesh points are in SI meters, and if not, warns the user, scales the
points to meters, and rewrites the mesh.

Intended to be called in `Allrun` scripts immediately after `newVtkUnstructuredToFoam`:

```bash
newVtkUnstructuredToFoam heart.vtk
checkMeshGeometry
cardiacFoam
```

## Behaviour

1. Read `constant/polyMesh` via `polyMesh`.
2. Compute bounding box and `maxDim = cmptMax(bb.span())`.
3. Apply auto-detection algorithm to determine `scaleFactor`.
4. If `scaleFactor == 1.0`: print summary, exit cleanly (no write).
5. If `scaleFactor != 1.0`: print `WarningInFunction` with detected units and
   scale factor, scale all points in-place, rewrite `constant/polyMesh`.
6. Always print: original bounding box, inferred unit, scale applied.

## Auto-detection algorithm

Target range for a cardiac mesh in meters: `maxDim ∈ [0.01, 1.0] m`.

```
scaleFactor = 10 ^ ( -1 - floor( log10(maxDim) ) )
```

| maxDim   | Inferred unit | scaleFactor |
|----------|---------------|-------------|
| ~0.1     | m (correct)   | 1 (no-op)   |
| ~100     | mm            | 1e-3        |
| ~350 000 | μm            | 1e-6        |

## File structure

```
applications/utilities/checkMeshGeometry/
├── checkMeshGeometry.C
└── Make/
    ├── files
    └── options
```

`Make/files`:

```
checkMeshGeometry.C
EXE = $(FOAM_USER_APPBIN)/checkMeshGeometry
```

`Make/options`: links `finiteVolume` and `meshTools`.

No separate header files — logic is self-contained in `checkMeshGeometry.C`.

`src/Allwmake` is updated to build this utility alongside the others.
