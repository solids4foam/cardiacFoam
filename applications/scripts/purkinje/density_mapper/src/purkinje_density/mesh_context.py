from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .vtk_bullseye_workbench import (
    SegmentFieldSelection,
    SegmentSurfaceAreaResult,
    compute_endocardial_surface_area_by_segment,
    discover_segment_field,
)


@dataclass(frozen=True)
class SegmentMeshContext:
    mesh: Any
    selection: SegmentFieldSelection
    segment_ids: Any
    area_result: SegmentSurfaceAreaResult


def extract_segment_ids(mesh: Any, selection: SegmentFieldSelection) -> Any:
    import numpy as np

    if selection.association == "point":
        return np.asarray(mesh.point_data[selection.field_name]).reshape(-1).astype(int)
    return np.asarray(mesh.cell_data[selection.field_name]).reshape(-1).astype(int)


def build_segment_mesh_context(
    vtk_path: Path,
    *,
    segment_field: str | None = None,
    segment_association: str = "auto",
    endocardial_selector_field: str | None = None,
    endocardial_selector_value: Any = 1,
) -> SegmentMeshContext:
    import pyvista as pv

    mesh = pv.read(str(vtk_path))
    selection = discover_segment_field(
        mesh=mesh,
        requested_field=segment_field,
        association=segment_association,
    )
    segment_ids = extract_segment_ids(mesh, selection)
    area_result = compute_endocardial_surface_area_by_segment(
        mesh=mesh,
        segment_selection=selection,
        segment_ids=segment_ids,
        endocardial_selector_field=endocardial_selector_field,
        endocardial_selector_value=endocardial_selector_value,
    )
    return SegmentMeshContext(
        mesh=mesh,
        selection=selection,
        segment_ids=segment_ids,
        area_result=area_result,
    )


def save_mesh_ascii(mesh: Any, output_path: str | Path) -> Path:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    mesh.save(str(out), binary=False)
    return out
