# newVtkUnstructuredToFoam

Extended VTK unstructured grid importer that reads scalar, vector, and tensor
CELL_DATA fields in addition to the mesh topology. Produces `constant/polyMesh`
and the imported fields in the current time directory (`0/` by default).

## Usage

```bash
newVtkUnstructuredToFoam myMesh.vtk -case ./myCase
```

Pass `-no-fields` to import mesh topology only without reading CELL_DATA.

## Complete post-import workflow

After import, two manual steps are always required before running cardiacFoam:

### Step 1 — Scale mesh to SI metres

VTK files from meshing tools (GMSH, Meshalyzer, SimNIBS, etc.) are typically in
millimetres. OpenFOAM requires SI metres.

```bash
transformPoints -scale '(0.001 0.001 0.001)' -case ./myCase
```

Verify with `checkMeshGeometry` — it reports the bounding box. A 20 × 3 × 7 cm
heart appears as ~0.2 x 0.03 x 0.07 m, not 200 x 30 x 70 mm.

### Step 2 — Fix conductivity field dimensions

All fields are written as dimensionless ([0 0 0 0 0 0 0]); the VTK format
carries no SI unit information. Most fields (fiber, sheet, tags,
uvc_transmural) are genuinely dimensionless and need no change.

conductivity is the conductivity tensor and must carry SI units S/m:

```bash
sed -i.bak 's/dimensions.*\[0 0 0 0 0 0 0\]/dimensions      [-1 -3 3 0 0 2 0]/' \
    ./myCase/0/conductivity
```

The correct dimension [-1 -3 3 0 0 2 0] encodes S/m (kg^-1 m^-3 s^3 A^2).

When to scale the values too: if the VTK file exported conductivity in S/mm
(not S/m), values must also be multiplied by 1000. Most cardiac modelling
pipelines export in S/m — check your meshing tool documentation.

## Dimension reference for common fields

| Field          | Class          | Correct dimensions   | Notes                          |
|----------------|----------------|----------------------|--------------------------------|
| conductivity   | volTensorField | [-1 -3 3 0 0 2 0]   | Conductivity tensor — fix this |
| fiber          | volVectorField | [0 0 0 0 0 0 0]     | Unit fibre direction — OK      |
| sheet          | volVectorField | [0 0 0 0 0 0 0]     | Unit sheet direction — OK      |
| uvc_transmural | volScalarField | [0 0 0 0 0 0 0]     | Transmural distance 0-1 — OK   |
| tags           | volScalarField | [0 0 0 0 0 0 0]     | Region labels — OK             |

## Quality checks after import

```bash
checkMeshGeometry -case ./myCase   # bounding box + unit sanity
checkMesh -case ./myCase           # topology, non-orthogonality, skewness
```

## Known limitations

- All CELL_DATA fields are written as dimless — see Step 2 above.
- No coordinate scaling — import is at source units; always run transformPoints.
- The script vtk_convert_arrays_to_fields.py referenced in older documentation
  does not exist in the repository.
