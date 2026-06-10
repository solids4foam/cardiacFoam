# setFieldDimensions

Assigns correct SI dimensions to field files that were produced by
`newVtkUnstructuredToFoam`, which writes every CELL_DATA field as
dimensionless `[0 0 0 0 0 0 0]`; the VTK format carries no unit information.

The utility patches **only the `dimensions` header line**. It never reads
or rewrites `internalField` or `boundaryField`, so it is safe and fast
even on very large meshes (millions of cells).

## Built-in catalogue

| Field                    | Applied dimensions    | Physical meaning      |
|--------------------------|-----------------------|-----------------------|
| `Diffusivity`            | `[-1 -3 3 0 0 2 0]`  | Conductivity tensor S/m |
| `conductivity`           | `[-1 -3 3 0 0 2 0]`  | Conductivity tensor S/m |
| `bodyAndOrgansConductivity` | `[-1 -3 3 0 0 2 0]` | Torso conductivity S/m |

Fields that are **not** in the catalogue (`fiber`, `sheet`, `uvc_transmural`,
`tags`) are genuinely dimensionless and are left unchanged.

## Usage

```bash
# Fix all catalogue fields in 0/ (typical post-VTK-import step)
setFieldDimensions -case ./myCase

# Fix a specific field
setFieldDimensions -field Diffusivity -case ./myCase

# Fix a custom field with explicit dimensions
setFieldDimensions -field myField -dim "[-1 -3 3 0 0 2 0]" -case ./myCase

# Process a different time directory
setFieldDimensions -time 0.1 -case ./myCase

# Preview without writing
setFieldDimensions -dryRun -case ./myCase
```

## Complete post-VTK-import sequence

```bash
newVtkUnstructuredToFoam  myHeart.vtk        -case ./myCase
transformPoints -scale '(0.001 0.001 0.001)' -case ./myCase
setFieldDimensions                           -case ./myCase
checkMeshGeometry                            -case ./myCase
checkMesh                                    -case ./myCase
```

## Build

```bash
wmake
```
