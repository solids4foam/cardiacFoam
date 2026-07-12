# checkMeshGeometry

Checks whether a mesh is in SI metres. Detection is read-only by default;
rescaling requires an explicit option.

The utility detects the unit by inspecting the maximum bounding-box dimension
of `constant/polyMesh`:

| Max dimension range | Detected unit | Automatic scale factor |
|---|---|---|
| \< 20 | m (already correct) | none |
| 20 – 999 | mm | 1e-3 |
| 1 000 – 999 999 | µm | 1e-6 |
| ≥ 1 000 000 | m (no automatic scaling) | none |

When rescaling is requested, the utility overwrites the selected region's
`points` and prints the original and scaled bounding boxes.

## Usage

```bash
checkMeshGeometry
checkMeshGeometry -rescale
checkMeshGeometry -scale 0.001
checkMeshGeometry -region <name>
```

## Options

- `-rescale` — apply the automatically detected factor and rewrite the mesh.
- `-scale <factor>` — apply an explicit factor and rewrite the mesh; this
  overrides automatic detection.
- `-region <name>` — operate on the named mesh region (default: `region0`).

## Notes

- Runs in serial only (`noParallel`).
- Without `-rescale` or `-scale`, the utility only reports its detection and
  does not modify the mesh.
- No other mesh files (boundary, faces, owner, neighbour) are modified; only
  `points` is rewritten via `polyMesh::movePoints` + `write()`.
- Run this utility before any solver that requires SI-unit coordinates.
