
# VTK Pre-processing and Field Conversion for OpenFOAM

This repository provides a Python-based workflow to process legacy ASCII VTK files for cardiac simulations, with a focus on generating OpenFOAM-compatible input files. It includes tools to add Purkinje layers, define diffusivity tensors, and convert generic VTK FIELD arrays into standard VTK data attributes.

Workflow Overview


 # vtk_pre_processing.py 
 Purkinje Layer and Diffusivity Tensor Generation with conversion to fields.

 The script reads a ASCIIlegacy.vtk and lets the user add a Purkinje network layer based on transmural distance and computes the diffusivity tensors for the cardiac mesh. Input for conductivity tensor and files written in <config.py>. It then runs <vtk_convert_arrays_to_fields.py> in the utils drectory to convert the input into OpenFOAM readable files. 

 This python script can be executed with:

```
python  vtk_pre_processing.py 
```

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
python  vtk_convert_arrays_to_fields.py
```


```
newVtkUnstructuredToFoam outputforFoam.vtk
```
 