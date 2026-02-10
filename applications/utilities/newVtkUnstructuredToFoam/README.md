# newVtkUnstructuredToFoam

`newVtkUnstructuredToFoam` is the same as `vtkUnstructuredToFoam` but has been
updated to also read and write tensor volFields.

It is recommended that you run the Python `vtk_convert_arrays_to_fields.py`
 script on the VTK file first to ensure all fields are in the correct format.
