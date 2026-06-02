# checkMeshGeometry

Checks whether a mesh is in SI meters and rescales it in-place if not.

The utility detects the unit by inspecting the maximum bounding-box dimension
of `constant/polyMesh`:

| Max dimension range | Detected unit | Scale factor applied |
|---|---|---|
| \< 1 | m (already correct) | none |
| 1 – 999 | mm | 1e-3 |
| 1 000 – 999 999 | µm | 1e-6 |

If rescaling is required, the utility overwrites `constant/polyMesh/points`
and prints the original and scaled bounding boxes.

## Usage

```bash
checkMeshGeometry
checkMeshGeometry -noScale    # detect and warn, but do not write
```

## Options

- `-noScale` — print the warning and scale factor but suppress the rewrite.

## Notes

- Runs in serial only (`noParallel`).
- The mesh is read from the default `constant/polyMesh` region.
- No other mesh files (boundary, faces, owner, neighbour) are modified; only
  `points` is rewritten via `polyMesh::movePoints` + `write()`.
- Run this utility before any solver that requires SI-unit coordinates.
