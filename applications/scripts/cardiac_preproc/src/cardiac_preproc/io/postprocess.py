"""Shared post-processing helpers for ASCII VTK outputs."""

from __future__ import annotations

from cardiac_preproc.utils.vtk_convert_arrays_to_fields import convert_vtk_file
from cardiac_preproc.utils.vtk_utils import inspect_fields, remove_blank_lines


def inspect_vtk(path: str, title: str | None = None) -> None:
    if title:
        print(f"\n-------------------------\n{title}\n-------------------------")
    inspect_fields(path)


def postprocess_vtk_output(
    path: str,
    *,
    convert_fields: bool = True,
    remove_blanks: bool = True,
) -> None:
    if not convert_fields and not remove_blanks:
        return
    print("\n-------------------------\nPost-processing VTK\n-------------------------")
    if remove_blanks:
        remove_blank_lines(path, path)
    if convert_fields:
        convert_vtk_file(path, path)
        print("Converted FIELD arrays to SCALARS/VECTORS/TENSORS.")
