"""Purkinje slab pipeline step."""

from __future__ import annotations

from dataclasses import dataclass

from cardiac_core.io.field_checks import require_point_fields
from cardiac_core.io.postprocess import inspect_vtk, postprocess_vtk_output
from cardiac_core.io.vtk_mesh import read_vtk_mesh, write_vtk_ascii
from cardiac_core.pipeline.context import StepContext, StepResult


@dataclass
class PurkinjeSlabOptions:
    input_path: str
    output_path: str
    transmural_min: float = 0.0
    transmural_max: float = 0.1
    lv_value: int = -1
    rv_value: int = 1
    field_name: str = "purkinjeLayer"
    purkinje_mult: float | None = None
    wall_tag_name: str | None = None
    wall_tag_value: float | int | None = None
    inflate_seed_point: str | None = None
    inflate_inside_tag_name: str = "inside_shared_boundary"
    inspect: bool = False
    convert_fields: bool = True
    remove_blank_lines: bool = True


def run_purkinje_slab(options: PurkinjeSlabOptions) -> StepResult:
    from cardiac_core.purkinje_network.purkinje_fractal.slab import (
        add_purkinje_layer,
    )

    if options.inspect:
        inspect_vtk(options.input_path, "Inspecting Input Fields")

    mesh = read_vtk_mesh(options.input_path)
    require_point_fields(
        mesh,
        ("uvc_transmural", "uvc_intraventricular"),
        "purkinje_slab",
    )

    mesh = add_purkinje_layer(
        mesh,
        transmural_min=options.transmural_min,
        transmural_max=options.transmural_max,
        lv_value=options.lv_value,
        rv_value=options.rv_value,
        field_name=options.field_name,
        purkinje_mult=options.purkinje_mult,
        wall_tag_name=options.wall_tag_name,
        wall_tag_value=options.wall_tag_value,
        inflate_seed_point=options.inflate_seed_point,
        inflate_inside_tag_name=options.inflate_inside_tag_name,
    )

    output_path = write_vtk_ascii(mesh, options.output_path)
    postprocess_vtk_output(
        output_path,
        convert_fields=options.convert_fields,
        remove_blanks=options.remove_blank_lines,
    )
    if options.inspect:
        inspect_vtk(output_path, "Inspecting Output Fields")
    print(f"Purkinje slab mesh written to {output_path}")
    return StepResult(output_mesh=output_path)


def run_from_context(context: StepContext) -> StepResult:
    opts = context.options
    options = PurkinjeSlabOptions(
        input_path=opts.get("input_path") or (context.current_mesh or ""),
        output_path=opts["output_path"],
        transmural_min=float(opts.get("transmural_min", 0.0)),
        transmural_max=float(opts.get("transmural_max", 0.1)),
        lv_value=int(opts.get("lv_value", -1)),
        rv_value=int(opts.get("rv_value", 1)),
        field_name=str(opts.get("field_name", "purkinjeLayer")),
        purkinje_mult=(
            None if opts.get("purkinje_mult") is None else float(opts.get("purkinje_mult"))
        ),
        wall_tag_name=opts.get("wall_tag_name"),
        wall_tag_value=opts.get("wall_tag_value"),
        inflate_seed_point=opts.get("inflate_seed_point"),
        inflate_inside_tag_name=str(
            opts.get("inflate_inside_tag_name", "inside_shared_boundary")
        ),
        inspect=bool(opts.get("inspect", False)),
        convert_fields=bool(opts.get("convert_fields", True)),
        remove_blank_lines=bool(opts.get("remove_blank_lines", True)),
    )
    return run_purkinje_slab(options)
