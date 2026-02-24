"""Purkinje slab tagging utilities.

Includes optional wall/boundary restriction using either:
- cell_data tags (e.g. endo_surface_lv == 1), or
- gmsh field_data + gmsh:physical cell_data mapping.
- seed-based growth (inflate from point until shared boundary).
"""

from __future__ import annotations

import numpy as np
import pyvista as pv

from cardiac_core.tagging.tag_endo_epi_surface import (
    map_surface_tags_to_volume,
    tag_inside_shared_boundary,
    tag_shared_boundary,
)


def _flat(values: np.ndarray) -> np.ndarray:
    return np.asarray(values).reshape(-1)


def _cell_ids_from_cell_tag(
    mesh: pv.DataSet,
    tag_name: str,
    tag_value: float | int | None = None,
) -> np.ndarray | None:
    if tag_name not in mesh.cell_data:
        return None
    values = _flat(mesh.cell_data[tag_name])
    if tag_value is None:
        mask = values != 0
    else:
        mask = np.isclose(values.astype(float), float(tag_value))
    return np.where(mask)[0]


def _cell_ids_from_field_data_tag(
    mesh: pv.DataSet,
    tag_name: str,
) -> np.ndarray | None:
    if tag_name not in mesh.field_data:
        return None
    if "gmsh:physical" not in mesh.cell_data:
        return None
    tag_entry = _flat(mesh.field_data[tag_name])
    if tag_entry.size == 0:
        return None
    physical_id = int(tag_entry[0])
    physical = _flat(mesh.cell_data["gmsh:physical"]).astype(int)
    return np.where(physical == physical_id)[0]


def _point_mask_from_cell_ids(mesh: pv.DataSet, cell_ids: np.ndarray) -> np.ndarray:
    mask = np.zeros(mesh.n_points, dtype=bool)
    for cid in cell_ids:
        ids = mesh.get_cell(int(cid)).point_ids
        mask[np.asarray(ids, dtype=int)] = True
    return mask


def _parse_seed_point(seed_point: str | tuple[float, float, float] | np.ndarray) -> np.ndarray:
    if isinstance(seed_point, np.ndarray):
        out = seed_point.astype(float).reshape(-1)
    elif isinstance(seed_point, tuple):
        out = np.asarray(seed_point, dtype=float).reshape(-1)
    else:
        parts = [p.strip() for p in str(seed_point).split(",") if p.strip()]
        out = np.asarray([float(p) for p in parts], dtype=float).reshape(-1)
    if out.size != 3:
        raise ValueError(f"Expected seed point with 3 values, got {out.size}: {seed_point}")
    return out


def resolve_wall_point_mask(
    mesh: pv.DataSet,
    wall_tag_name: str,
    wall_tag_value: float | int | None = None,
) -> np.ndarray:
    """Return a boolean point mask for the selected wall tag."""
    cell_ids = _cell_ids_from_cell_tag(mesh, wall_tag_name, wall_tag_value)
    if cell_ids is None:
        cell_ids = _cell_ids_from_field_data_tag(mesh, wall_tag_name)
    if cell_ids is None:
        raise RuntimeError(
            "Unable to resolve wall tag "
            f"'{wall_tag_name}'. Expected it in cell_data or in field_data "
            "with companion 'gmsh:physical' cell_data."
        )
    if cell_ids.size == 0:
        raise RuntimeError(f"Wall tag '{wall_tag_name}' resolved to zero cells.")
    return _point_mask_from_cell_ids(mesh, cell_ids)


def resolve_inflate_point_mask(
    mesh: pv.DataSet,
    seed_point: str | tuple[float, float, float] | np.ndarray,
    *,
    intraventricular_field: str = "uvc_intraventricular",
    transmural_field: str = "uvc_transmural",
    lv_value: int = -1,
    rv_value: int = 1,
    rv_transmural_eq: float = 0.0,
    rv_transmural_eps: float = 1e-6,
    lv_transmural_eq: float = 1.0,
    lv_transmural_eps: float = 1e-6,
    shared_tag_name: str = "shared_boundary",
    inside_tag_name: str = "inside_shared_boundary",
) -> np.ndarray:
    """Grow from a seed on LV epi until shared boundary and return a point mask."""
    seed = _parse_seed_point(seed_point)
    surface = mesh.extract_surface().triangulate()
    surface = tag_shared_boundary(
        surface=surface,
        transmural_field=transmural_field,
        intraventricular_field=intraventricular_field,
        rv_transmural_eq=rv_transmural_eq,
        rv_transmural_eps=rv_transmural_eps,
        lv_transmural_eq=lv_transmural_eq,
        lv_transmural_eps=lv_transmural_eps,
        lv_value=lv_value,
        rv_value=rv_value,
        tag_name=shared_tag_name,
    )
    surface = tag_inside_shared_boundary(
        surface=surface,
        seed_point=seed,
        intraventricular_field=intraventricular_field,
        transmural_field=transmural_field,
        lv_value=lv_value,
        lv_transmural_eq=lv_transmural_eq,
        lv_transmural_eps=lv_transmural_eps,
        shared_tag_name=shared_tag_name,
        inside_tag_name=inside_tag_name,
    )

    volume = map_surface_tags_to_volume(mesh.copy(deep=True), surface, [inside_tag_name])
    if inside_tag_name not in volume.cell_data:
        raise RuntimeError(
            f"Failed to map '{inside_tag_name}' from surface growth back to volume."
        )
    inside_cells = np.where(_flat(volume.cell_data[inside_tag_name]) > 0)[0]
    if inside_cells.size == 0:
        raise RuntimeError(
            f"Inflate-from-point produced zero cells for tag '{inside_tag_name}'."
        )
    return _point_mask_from_cell_ids(mesh, inside_cells)


def add_purkinje_layer(
    mesh: pv.DataSet,
    transmural_min: float = 0.0,
    transmural_max: float = 0.1,
    lv_value: int = -1,
    rv_value: int = 1,
    field_name: str = "purkinjeLayer",
    purkinje_mult: float | None = None,
    diffusivity_field: str = "Diffusivity",
    wall_tag_name: str | None = None,
    wall_tag_value: float | int | None = None,
    inflate_seed_point: str | tuple[float, float, float] | np.ndarray | None = None,
    inflate_inside_tag_name: str = "inside_shared_boundary",
) -> pv.DataSet:
    """Tag a Purkinje slab from UVC fields with optional wall restriction."""
    if "uvc_transmural" not in mesh.point_data or "uvc_intraventricular" not in mesh.point_data:
        raise RuntimeError(
            "Missing required point_data fields 'uvc_transmural' and 'uvc_intraventricular'."
        )

    uvc_transmural = _flat(mesh.point_data["uvc_transmural"])
    uvc_intraventricular = _flat(mesh.point_data["uvc_intraventricular"])

    is_endo = (uvc_transmural >= transmural_min) & (uvc_transmural <= transmural_max)
    is_rv = uvc_intraventricular == rv_value
    is_lv = uvc_intraventricular == lv_value
    purkinje_labels = is_endo & (is_rv | is_lv)

    if inflate_seed_point is not None:
        inflate_mask = resolve_inflate_point_mask(
            mesh,
            inflate_seed_point,
            lv_value=lv_value,
            rv_value=rv_value,
            inside_tag_name=inflate_inside_tag_name,
        )
        purkinje_labels |= inflate_mask
        print(
            "Applied inflate-from-point additive region "
            f"(points added candidate-set: {int(np.count_nonzero(inflate_mask))})."
        )

    if wall_tag_name:
        wall_point_mask = resolve_wall_point_mask(mesh, wall_tag_name, wall_tag_value)
        purkinje_labels &= wall_point_mask
        print(
            f"Applied wall filter '{wall_tag_name}' "
            f"(points kept: {int(np.count_nonzero(wall_point_mask))})."
        )

    cell_tags = np.zeros(mesh.n_cells, dtype=np.int8)
    for i in range(mesh.n_cells):
        point_ids = mesh.get_cell(i).point_ids
        if np.any(purkinje_labels[np.asarray(point_ids, dtype=int)]):
            cell_tags[i] = 1

    mesh.cell_data[field_name] = cell_tags
    percent = (transmural_max - transmural_min) * 100.0
    print(
        f"Added purkinje layer in {percent:g}% of transmural depth "
        f"({np.count_nonzero(cell_tags)} cells tagged as {field_name})."
    )

    if purkinje_mult is not None and diffusivity_field in mesh.cell_data:
        diffusivity = np.asarray(mesh.cell_data[diffusivity_field])
        if diffusivity.ndim != 2 or diffusivity.shape[1] != 9:
            raise RuntimeError(
                f"Expected {diffusivity_field} to be (n_cells, 9), got {diffusivity.shape}."
            )
        mask = cell_tags == 1
        diffusivity = diffusivity.copy()
        diffusivity[mask] *= purkinje_mult
        mesh.cell_data[diffusivity_field] = diffusivity
        print(
            f"Scaled {diffusivity_field} by {purkinje_mult}x on "
            f"{np.count_nonzero(mask)} purkinjeLayer cells."
        )
    elif purkinje_mult is not None:
        print(
            f"Warning: {diffusivity_field} not found; purkinjeLayer tagged without scaling."
        )

    return mesh
