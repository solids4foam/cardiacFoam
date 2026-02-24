"""Shared VTK mesh read/write helpers for preprocessing steps."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def read_vtk_mesh(path: str) -> Any:
    import pyvista as pv

    return pv.read(path)


def write_vtk_ascii(mesh: Any, path: str) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    mesh.save(str(output), binary=False)
    return str(output)
