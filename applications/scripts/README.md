# vtk_convert_arrays_to_fields.py

This Python script processes legacy ASCII VTK files and converts generic FIELD
 arrays into standard VTK data attributes:

- 1-component fields → SCALARS <name> <type> (with LOOKUP_TABLE default)
- 3-component fields → VECTORS <name> <type>
- 9-component fields → TENSORS <name> <type> (written as 3×3 blocks)

The script scans `POINT_DATA` and `CELL_DATA` sections, replaces convertible
 arrays with their appropriate VTK keyword, and rebuilds or removes the `FIELD`
 blocks as needed.

This ensures scalar, vector, and tensor data are recognised natively by VTK/ParaView
 without needing manual selection from `FieldData`. Also, it allows the fields to
 be correctly parsed by the OpenFOAM `newVtkUnstructuredToFoam` conversion
 utility.

This python script can be executed with:

```
python3 vtk_convert_arrays_to_fields.py
```

The script expects to find a `ASCIIlegacy.vtk` file in the current directory,
 and will write a new file called `ASCIIlegacy_fixed.vtk`. The `ASCIIlegacy_fixed.vtk`
 file can be converted to OpenFOAM mesh and fields with:

```
newVtkUnstructuredToFoam ASCIIlegacy_fixed.vtk
```
 