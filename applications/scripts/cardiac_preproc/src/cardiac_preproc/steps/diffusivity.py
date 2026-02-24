"""Diffusivity pipeline step."""

from __future__ import annotations

from dataclasses import dataclass

from cardiac_preproc.io.field_checks import require_cell_fields
from cardiac_preproc.io.vtk_mesh import read_vtk_mesh, write_vtk_ascii
from cardiac_preproc.io.postprocess import inspect_vtk, postprocess_vtk_output
from cardiac_preproc.lib.diffusivity_tensor import add_diffusivity_tensor_ventricles
from cardiac_preproc.pipeline.context import StepContext, StepResult


@dataclass
class DiffusivityOptions:
    input_path: str
    output_path: str
    df: float
    ds: float
    dn: float
    purkinje_mult: float | None = None
    inspect: bool = True
    convert_fields: bool = True
    remove_blank_lines: bool = True


def run_diffusivity(options: DiffusivityOptions) -> StepResult:
    if options.inspect:
        inspect_vtk(options.input_path, "Inspecting Input Fields")

    mesh = read_vtk_mesh(options.input_path)
    require_cell_fields(mesh, ("fiber", "sheet"), "diffusivity")
    has_purkinje = "purkinjeLayer" in mesh.cell_data
    if has_purkinje:
        print(
            "Warning: purkinjeLayer found in input. "
            "Preferred order is diffusivity -> purkinje_slab."
        )

    scale = 1.0 if options.purkinje_mult is None else options.purkinje_mult
    mesh = add_diffusivity_tensor_ventricles(
        mesh,
        options.df,
        options.ds,
        options.dn,
        purkinje_mult=scale,
    )
    write_vtk_ascii(mesh, options.output_path)

    postprocess_vtk_output(
        options.output_path,
        convert_fields=options.convert_fields,
        remove_blanks=options.remove_blank_lines,
    )
    if options.inspect:
        inspect_vtk(options.output_path, "Inspecting Output Fields")
    return StepResult(output_mesh=options.output_path)


def run_from_context(context: StepContext) -> StepResult:
    opts = context.options
    options = DiffusivityOptions(
        input_path=opts.get("input_path") or (context.current_mesh or ""),
        output_path=opts["output_path"],
        df=float(opts["df"]),
        ds=float(opts["ds"]),
        dn=float(opts["dn"]),
        purkinje_mult=opts.get("purkinje_mult"),
        inspect=bool(opts.get("inspect", False)),
        convert_fields=bool(opts.get("convert_fields", True)),
        remove_blank_lines=bool(opts.get("remove_blank_lines", True)),
    )
    return run_diffusivity(options)
