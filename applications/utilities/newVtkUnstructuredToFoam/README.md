# newVtkUnstructuredToFoam

Extended VTK unstructured grid importer that reads scalar, vector, and tensor
CELL_DATA fields in addition to the mesh topology. Produces `constant/polyMesh`
and the imported fields in the current time directory (`0/` by default).

## Usage

```bash
newVtkUnstructuredToFoam myMesh.vtk -case ./myCase
```

## After import

Imported coordinates are at source units, and every CELL_DATA field is
written dimensionless — VTK carries no unit information. Scaling the mesh to
SI metres and fixing field dimensions are always required next; see
[utilities: Importing an external mesh](../README.md#importing-an-external-mesh)
for the full sequence and which tool owns each step.

## Known limitations

- Import supports the documented VTK subset only; unsupported CELL_DATA
  types are skipped.
- No coordinate scaling: import is at source units.
