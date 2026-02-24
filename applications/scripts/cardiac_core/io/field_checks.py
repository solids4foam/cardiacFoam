"""Reusable mesh field checks for pipeline steps."""

from __future__ import annotations

from collections.abc import Iterable

def require_cell_fields(mesh, required: Iterable[str], step_name: str) -> None:
    missing = [name for name in required if name not in mesh.cell_data]
    if missing:
        raise RuntimeError(
            f"{step_name}: missing required CELL_DATA field(s): {', '.join(missing)}"
        )


def require_point_fields(mesh, required: Iterable[str], step_name: str) -> None:
    missing = [name for name in required if name not in mesh.point_data]
    if missing:
        raise RuntimeError(
            f"{step_name}: missing required POINT_DATA field(s): {', '.join(missing)}"
        )
