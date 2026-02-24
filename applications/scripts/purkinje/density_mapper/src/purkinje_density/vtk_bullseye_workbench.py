from __future__ import annotations

import argparse
import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence
import sys

try:
    import numpy as np  # type: ignore
except Exception:
    np = None  # type: ignore

if __package__ in {None, ""}:
    # Support direct script execution:
    #   python src/purkinje_density/vtk_bullseye_workbench.py ...
    # by adding `src/` to sys.path and using absolute imports.
    _src_root = Path(__file__).resolve().parents[1]
    if str(_src_root) not in sys.path:
        sys.path.insert(0, str(_src_root))
    from purkinje_density.interactive_bullseye import (  # type: ignore
        DEFAULT_AHA17_REFERENCE_CFG,
        collect_segment_percentages_interactive,
        segment_ids_from_reference,
    )
else:
    from .interactive_bullseye import (
        DEFAULT_AHA17_REFERENCE_CFG,
        collect_segment_percentages_interactive,
        segment_ids_from_reference,
    )

DEFAULT_SEGMENT_FIELD_CANDIDATES: tuple[str, ...] = (
    "anatomical_tag",
    "tagged_segment_id",
    "anatomical_segment_id",
    "segment_id",
    "seg_id",
)
DEFAULT_ENDOCARDIAL_SELECTOR_FIELD_CANDIDATES: tuple[str, ...] = (
    "endo_surface",
    "endo_surface_lv",
    "endo_surface_rv",
)


@dataclass(frozen=True)
class SegmentFieldSelection:
    field_name: str
    association: str  # "point" or "cell"


@dataclass(frozen=True)
class DivisionReference:
    reference_name: str
    reference_cfg: dict[str, object]
    source: str


@dataclass(frozen=True)
class SegmentSurfaceAreaResult:
    area_by_segment: dict[int, float]
    mapped_values: Any
    selector_source: str


def _load_python_module(module_path: Path, module_name: str) -> Any:
    spec = importlib.util.spec_from_file_location(module_name, str(module_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module from '{module_path}'.")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _default_anatomical_config_path() -> Path:
    # .../purkinje/density_mapper/src/purkinje_density/vtk_bullseye_workbench.py
    # -> .../applications/scripts/anatomicalSegmentationModels/configs/lv_division_config.py
    return (
        Path(__file__).resolve().parents[3]
        / "anatomicalSegmentationModels"
        / "configs"
        / "lv_division_config.py"
    )


def load_division_reference(
    reference_name: str | None = None,
    config_source: str = "anatomical",
    config_file: Path | None = None,
) -> DivisionReference:
    """
    Load bullseye division configuration.

    `config_source`:
    - anatomical: from anatomicalSegmentationModels config file.
    - local: built-in AHA17 fallback.
    """
    source = config_source.strip().lower()
    if source not in {"anatomical", "local"}:
        raise ValueError("config_source must be 'anatomical' or 'local'.")

    if source == "local":
        ref_name = reference_name or "aha17"
        if ref_name != "aha17":
            raise KeyError(
                f"Local config only provides 'aha17', got '{ref_name}'."
            )
        return DivisionReference(
            reference_name="aha17",
            reference_cfg=dict(DEFAULT_AHA17_REFERENCE_CFG),
            source="local",
        )

    cfg_path = config_file or _default_anatomical_config_path()
    if not cfg_path.exists():
        if reference_name is not None and reference_name != "aha17":
            raise FileNotFoundError(
                f"Anatomical config file '{cfg_path}' not found and reference "
                f"'{reference_name}' is not available locally."
            )
        return DivisionReference(
            reference_name="aha17",
            reference_cfg=dict(DEFAULT_AHA17_REFERENCE_CFG),
            source="local-fallback",
        )

    mod = _load_python_module(cfg_path, "lv_division_config_ext")
    references = dict(getattr(mod, "division_references", {}))
    active_ref = reference_name or str(getattr(mod, "division_reference", "")).strip()
    if active_ref == "":
        active_ref = "aha17"

    if active_ref not in references:
        available = ", ".join(sorted(references.keys()))
        raise KeyError(
            f"Division reference '{active_ref}' not found in '{cfg_path}'. "
            f"Available: {available}"
        )

    ref_cfg = dict(references[active_ref])
    return DivisionReference(
        reference_name=str(active_ref),
        reference_cfg=ref_cfg,
        source=str(cfg_path),
    )


def map_percentages_to_segments(
    segment_ids: Sequence[int] | Any,
    segment_percentages: Mapping[int, float],
    default_value: float = 0.0,
    ignore_non_positive: bool = True,
) -> Any:
    """
    Map segment IDs to percentage values.
    """
    if np is None:
        ids = [int(x) for x in segment_ids]
        out: list[float] = []
        for seg_id in ids:
            if ignore_non_positive and int(seg_id) <= 0:
                out.append(float(default_value))
            else:
                out.append(float(segment_percentages.get(int(seg_id), default_value)))
        return out

    ids = np.asarray(segment_ids).reshape(-1).astype(int)
    mapped = np.full(ids.shape, float(default_value), dtype=float)
    if mapped.size == 0:
        return mapped

    unique_ids = np.unique(ids)
    for seg_id in unique_ids:
        if ignore_non_positive and int(seg_id) <= 0:
            continue
        value = float(segment_percentages.get(int(seg_id), default_value))
        mapped[ids == int(seg_id)] = value
    return mapped


def _selector_mask(values: Any, selector_value: Any) -> Any:
    if np is None:
        return [str(v) == str(selector_value) for v in values]

    arr = np.asarray(values).reshape(-1)
    if arr.size == 0:
        return np.zeros((0,), dtype=bool)

    if arr.dtype.kind in {"U", "S", "O"}:
        return np.asarray([str(v) == str(selector_value) for v in arr], dtype=bool)

    try:
        target = float(selector_value)
        return np.isclose(arr.astype(float), target)
    except Exception:
        return np.asarray([str(v) == str(selector_value) for v in arr], dtype=bool)


def aggregate_surface_area_by_segment(
    cell_segment_ids: Sequence[int] | Any,
    cell_areas: Sequence[float] | Any,
    cell_mask: Sequence[bool] | Any | None = None,
    ignore_non_positive: bool = True,
) -> dict[int, float]:
    if np is None:
        out: dict[int, float] = {}
        if cell_mask is None:
            cell_mask = [True] * len(cell_segment_ids)
        for seg_id, area, keep in zip(cell_segment_ids, cell_areas, cell_mask):
            seg = int(seg_id)
            if not keep:
                continue
            if ignore_non_positive and seg <= 0:
                continue
            out[seg] = float(out.get(seg, 0.0) + float(area))
        return out

    seg = np.asarray(cell_segment_ids).reshape(-1).astype(int)
    areas = np.asarray(cell_areas).reshape(-1).astype(float)
    if seg.shape[0] != areas.shape[0]:
        raise ValueError("cell_segment_ids and cell_areas must have the same length.")

    if cell_mask is None:
        mask = np.ones(seg.shape[0], dtype=bool)
    else:
        mask = np.asarray(cell_mask).reshape(-1).astype(bool)
        if mask.shape[0] != seg.shape[0]:
            raise ValueError("cell_mask must have the same length as cell_segment_ids.")

    valid = mask
    if ignore_non_positive:
        valid = valid & (seg > 0)

    out: dict[int, float] = {}
    for seg_id in np.unique(seg[valid]):
        out[int(seg_id)] = float(areas[valid & (seg == int(seg_id))].sum())
    return out


def _cell_segment_ids_from_point_segments(mesh: Any, point_segment_ids: Any) -> Any:
    if np is None:
        out: list[int] = []
        for cell_id in range(mesh.n_cells):
            pids = list(mesh.get_cell(cell_id).point_ids)
            vals = [int(point_segment_ids[i]) for i in pids if int(point_segment_ids[i]) > 0]
            if not vals:
                out.append(0)
                continue
            freq: dict[int, int] = {}
            for v in vals:
                freq[v] = freq.get(v, 0) + 1
            out.append(max(freq.items(), key=lambda kv: kv[1])[0])
        return out

    point_ids = np.asarray(point_segment_ids).reshape(-1).astype(int)
    out = np.zeros((mesh.n_cells,), dtype=int)
    for cell_id in range(mesh.n_cells):
        pids = np.asarray(mesh.get_cell(cell_id).point_ids).reshape(-1).astype(int)
        vals = point_ids[pids]
        vals = vals[vals > 0]
        if vals.size == 0:
            out[cell_id] = 0
            continue
        uniq, cnt = np.unique(vals, return_counts=True)
        out[cell_id] = int(uniq[np.argmax(cnt)])
    return out


def _cell_mask_from_point_mask(
    mesh: Any,
    point_mask: Any,
    threshold: float = 0.5,
) -> Any:
    if np is None:
        out: list[bool] = []
        for cell_id in range(mesh.n_cells):
            pids = list(mesh.get_cell(cell_id).point_ids)
            if not pids:
                out.append(False)
                continue
            frac = sum(1 for i in pids if point_mask[i]) / float(len(pids))
            out.append(frac >= float(threshold))
        return out

    pmask = np.asarray(point_mask).reshape(-1).astype(bool)
    out = np.zeros((mesh.n_cells,), dtype=bool)
    for cell_id in range(mesh.n_cells):
        pids = np.asarray(mesh.get_cell(cell_id).point_ids).reshape(-1).astype(int)
        if pids.size == 0:
            out[cell_id] = False
            continue
        out[cell_id] = bool(pmask[pids].mean() >= float(threshold))
    return out


def _mask_for_selector_field(
    mesh: Any,
    field_name: str,
    selector_value: Any,
) -> Any:
    if field_name in mesh.cell_data:
        return _selector_mask(mesh.cell_data[field_name], selector_value)
    if field_name in mesh.point_data:
        point_mask = _selector_mask(mesh.point_data[field_name], selector_value)
        return _cell_mask_from_point_mask(mesh, point_mask, threshold=0.5)
    raise KeyError(f"Selector field '{field_name}' not found in point_data or cell_data.")


def _auto_endocardial_cell_mask(
    mesh: Any,
    selector_value: Any,
) -> tuple[Any, str] | None:
    if np is None:
        return None

    primary = DEFAULT_ENDOCARDIAL_SELECTOR_FIELD_CANDIDATES[0]
    if primary in mesh.cell_data or primary in mesh.point_data:
        return (
            _mask_for_selector_field(mesh, primary, selector_value),
            f"auto:{primary}",
        )

    subfields = [
        f
        for f in DEFAULT_ENDOCARDIAL_SELECTOR_FIELD_CANDIDATES[1:]
        if (f in mesh.cell_data or f in mesh.point_data)
    ]
    if not subfields:
        return None

    combined = np.zeros((mesh.n_cells,), dtype=bool)
    for f in subfields:
        combined = combined | np.asarray(_mask_for_selector_field(mesh, f, selector_value)).reshape(-1).astype(bool)
    return combined, f"auto:{'+'.join(subfields)}"


def compute_endocardial_surface_area_by_segment(
    mesh: Any,
    segment_selection: SegmentFieldSelection,
    segment_ids: Sequence[int] | Any,
    endocardial_selector_field: str | None = None,
    endocardial_selector_value: Any = 1,
) -> SegmentSurfaceAreaResult:
    """
    Compute endocardial surface area per segment and map it back to current
    segment association (point/cell) as a scalar array.

    If `endocardial_selector_field` is None, tries auto-detection using
    DEFAULT_ENDOCARDIAL_SELECTOR_FIELD_CANDIDATES. If none are present, all
    mesh cells are treated as endocardial surface cells.
    """
    if np is None:
        raise RuntimeError("NumPy is required to compute surface areas.")

    sized = mesh.compute_cell_sizes(length=False, area=True, volume=False)
    if "Area" not in sized.cell_data:
        raise KeyError("Could not compute cell area field 'Area'.")
    cell_areas = np.asarray(sized.cell_data["Area"]).reshape(-1).astype(float)

    if segment_selection.association == "cell":
        cell_segment_ids = np.asarray(mesh.cell_data[segment_selection.field_name]).reshape(-1).astype(int)
    else:
        point_segment_ids = np.asarray(mesh.point_data[segment_selection.field_name]).reshape(-1).astype(int)
        cell_segment_ids = _cell_segment_ids_from_point_segments(mesh, point_segment_ids)

    if endocardial_selector_field is None:
        auto_result = _auto_endocardial_cell_mask(mesh, endocardial_selector_value)
        if auto_result is None:
            cell_mask = np.ones((mesh.n_cells,), dtype=bool)
            selector_source = "fallback:all_cells"
        else:
            cell_mask, selector_source = auto_result
    elif endocardial_selector_field in mesh.cell_data or endocardial_selector_field in mesh.point_data:
        cell_mask = _mask_for_selector_field(mesh, endocardial_selector_field, endocardial_selector_value)
        selector_source = f"explicit:{endocardial_selector_field}"
    else:
        raise KeyError(
            f"Endocardial selector field '{endocardial_selector_field}' not found "
            f"in point_data or cell_data."
        )

    area_by_segment = aggregate_surface_area_by_segment(
        cell_segment_ids=cell_segment_ids,
        cell_areas=cell_areas,
        cell_mask=cell_mask,
        ignore_non_positive=True,
    )
    mapped = map_percentages_to_segments(
        segment_ids=segment_ids,
        segment_percentages=area_by_segment,
        default_value=0.0,
    )
    return SegmentSurfaceAreaResult(
        area_by_segment=area_by_segment,
        mapped_values=mapped,
        selector_source=selector_source,
    )


def discover_segment_field(
    mesh: Any,
    requested_field: str | None = None,
    association: str = "auto",
) -> SegmentFieldSelection:
    association = association.strip().lower()
    if association not in {"auto", "point", "cell"}:
        raise ValueError("association must be one of: auto, point, cell")

    point_fields = set(mesh.point_data.keys())
    cell_fields = set(mesh.cell_data.keys())

    candidates: list[str] = []
    if requested_field:
        candidates.append(requested_field)
    for name in DEFAULT_SEGMENT_FIELD_CANDIDATES:
        if name not in candidates:
            candidates.append(name)

    if association in {"auto", "point"}:
        for name in candidates:
            if name in point_fields:
                return SegmentFieldSelection(field_name=name, association="point")
    if association in {"auto", "cell"}:
        for name in candidates:
            if name in cell_fields:
                return SegmentFieldSelection(field_name=name, association="cell")

    raise KeyError(
        "Could not find anatomical segment field in VTK. "
        f"Requested='{requested_field}', point fields={sorted(point_fields)}, "
        f"cell fields={sorted(cell_fields)}. "
        "If needed, run anatomical tagging first (anatomicalSegmentationModels) "
        "and then re-run this workbench."
    )


def _show_mesh(
    mesh: Any,
    scalars: str,
    preference: str,
    title: str,
    cmap: str,
    clim: tuple[float, float] | None = None,
) -> None:
    try:
        import pyvista as pv
    except Exception as exc:
        raise RuntimeError("PyVista is required for 3D visualization.") from exc

    pl = pv.Plotter()
    pl.add_mesh(
        mesh,
        scalars=scalars,
        preference=preference,
        cmap=cmap,
        clim=clim,
        show_edges=False,
    )
    pl.add_text(title, position="upper_left", font_size=11)
    pl.show()


def run_vtk_bullseye_workbench(
    vtk_path: Path,
    output_vtk_path: Path | None = None,
    reference_name: str | None = None,
    config_source: str = "anatomical",
    config_file: Path | None = None,
    segment_field: str | None = None,
    segment_association: str = "auto",
    mapped_field_name: str = "bullseye_percentage",
    show_original_3d: bool = True,
    show_mapped_3d: bool = True,
    prompt_missing_after_close: bool = True,
    input_mode: str = "widget",
    value_max: float = 10.0,
) -> dict[int, float]:
    try:
        import pyvista as pv
    except Exception as exc:
        raise RuntimeError("PyVista is required for VTK workbench.") from exc

    if not vtk_path.exists():
        raise FileNotFoundError(f"VTK file not found: {vtk_path}")

    division = load_division_reference(
        reference_name=reference_name,
        config_source=config_source,
        config_file=config_file,
    )

    mesh = pv.read(str(vtk_path))
    selection = discover_segment_field(
        mesh=mesh,
        requested_field=segment_field,
        association=segment_association,
    )
    if selection.association == "point":
        segment_ids_raw = mesh.point_data[selection.field_name]
    else:
        segment_ids_raw = mesh.cell_data[selection.field_name]

    if np is not None:
        segment_ids = np.asarray(segment_ids_raw).reshape(-1).astype(int)
        present_ids = sorted(int(x) for x in np.unique(segment_ids) if int(x) > 0)
    else:
        segment_ids = [int(x) for x in segment_ids_raw]
        present_ids = sorted({int(x) for x in segment_ids if int(x) > 0})

    expected_ids = sorted(segment_ids_from_reference(division.reference_cfg))
    unknown_ids = sorted(set(present_ids) - set(expected_ids))
    if unknown_ids:
        print(
            "Warning: VTK contains segment IDs not in selected bullseye reference: "
            f"{unknown_ids}"
        )

    if show_original_3d:
        _show_mesh(
            mesh=mesh,
            scalars=selection.field_name,
            preference=selection.association,
            title=(
                f"Original Anatomical Segments\n"
                f"field={selection.field_name} ({selection.association})"
            ),
            cmap="tab20",
        )

    percentages = collect_segment_percentages_interactive(
        reference_cfg=division.reference_cfg,
        title=f"Bullseye Input ({division.reference_name})",
        prompt_missing_after_close=prompt_missing_after_close,
        input_mode=input_mode,
        max_percentage_value=value_max,
    )

    mapped = map_percentages_to_segments(segment_ids, percentages, default_value=0.0)
    out_mesh = mesh.copy(deep=True)
    if selection.association == "point":
        out_mesh.point_data[mapped_field_name] = mapped
    else:
        out_mesh.cell_data[mapped_field_name] = mapped

    if output_vtk_path is not None:
        output_vtk_path.parent.mkdir(parents=True, exist_ok=True)
        out_mesh.save(str(output_vtk_path), binary=False)
        print(f"Saved mapped VTK: {output_vtk_path}")

    if show_mapped_3d:
        _show_mesh(
            mesh=out_mesh,
            scalars=mapped_field_name,
            preference=selection.association,
            title=(
                f"Mapped Bullseye Percentages\n"
                f"field={mapped_field_name} ({selection.association})"
            ),
            cmap="YlOrRd",
            clim=(0.0, float(value_max)),
        )

    return percentages


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Open VTK + interactive bullseye input, map segment percentages "
            "back to the mesh, and optionally save."
        )
    )
    ap.add_argument("--vtk", type=Path, required=True, help="Input VTK/VTU/VTP mesh path.")
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output VTK path with mapped percentage field.",
    )
    ap.add_argument(
        "--reference",
        default=None,
        help="Division reference name (for example: aha17, apical4_mid5_basal5).",
    )
    ap.add_argument(
        "--config-source",
        choices=["anatomical", "local"],
        default="anatomical",
        help="Where to load bullseye config from.",
    )
    ap.add_argument(
        "--config-file",
        type=Path,
        default=None,
        help="Optional path to a lv_division_config.py file.",
    )
    ap.add_argument(
        "--segment-field",
        default=None,
        help=(
            "Anatomical segment field name in VTK. If omitted, auto-detected from "
            f"{DEFAULT_SEGMENT_FIELD_CANDIDATES}."
        ),
    )
    ap.add_argument(
        "--segment-association",
        choices=["auto", "point", "cell"],
        default="auto",
        help="Whether segment IDs are in point_data or cell_data.",
    )
    ap.add_argument(
        "--mapped-field",
        default="bullseye_percentage",
        help="Output scalar field name for mapped percentages.",
    )
    ap.add_argument(
        "--no-show-original",
        action="store_true",
        help="Skip opening the original 3D mesh window.",
    )
    ap.add_argument(
        "--no-show-mapped",
        action="store_true",
        help="Skip opening the mapped 3D mesh window.",
    )
    ap.add_argument(
        "--allow-partial",
        action="store_true",
        help=(
            "Do not prompt for missing segments after bullseye closes; "
            "return only provided entries."
        ),
    )
    ap.add_argument(
        "--input-mode",
        choices=["widget", "terminal"],
        default="widget",
        help="Bullseye value-entry mode (widget is in-figure panel input).",
    )
    ap.add_argument(
        "--value-max",
        type=float,
        default=10.0,
        help="Maximum per-segment value for entry/color scale (default: 10).",
    )
    return ap


def main() -> None:
    ap = build_arg_parser()
    args = ap.parse_args()
    percentages = run_vtk_bullseye_workbench(
        vtk_path=args.vtk,
        output_vtk_path=args.output,
        reference_name=args.reference,
        config_source=args.config_source,
        config_file=args.config_file,
        segment_field=args.segment_field,
        segment_association=args.segment_association,
        mapped_field_name=args.mapped_field,
        show_original_3d=not args.no_show_original,
        show_mapped_3d=not args.no_show_mapped,
        prompt_missing_after_close=not args.allow_partial,
        input_mode=args.input_mode,
        value_max=args.value_max,
    )
    print("Final segment percentage map:")
    for seg_id in sorted(percentages.keys()):
        print(f"  {seg_id}: {percentages[seg_id]:.1f}%")


if __name__ == "__main__":
    main()
