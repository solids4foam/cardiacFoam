"""Scar pipeline step."""

from __future__ import annotations

from dataclasses import dataclass

from cardiac_preproc.io.vtk_mesh import read_vtk_mesh, write_vtk_ascii
from cardiac_preproc.io.postprocess import inspect_vtk, postprocess_vtk_output
from cardiac_preproc.lib.scar import add_scar_from_selection
from cardiac_preproc.pipeline.context import StepContext, StepResult


@dataclass
class ScarOptions:
    full_mesh_path: str
    selection_path: str
    output_path: str
    id_array: str = "GlobalCellIds"
    scar_name: str = "Scar"
    scar_value: float = 1.0
    diffusivity_name: str = "Diffusivity"
    diffusivity_scale: float = 0.0
    inspect: bool = False
    convert_fields: bool = True
    remove_blank_lines: bool = True


def run_scar(options: ScarOptions) -> StepResult:
    if options.inspect:
        inspect_vtk(options.full_mesh_path, "Inspecting Full Mesh Fields")
        inspect_vtk(options.selection_path, "Inspecting Selection Fields")

    full = read_vtk_mesh(options.full_mesh_path)
    sel = read_vtk_mesh(options.selection_path)

    full, scar_mask = add_scar_from_selection(
        full,
        sel,
        id_array=options.id_array,
        scar_name=options.scar_name,
        diffusivity_name=options.diffusivity_name,
        scar_value=options.scar_value,
        diffusivity_scale=options.diffusivity_scale,
        return_mask=True,
    )
    print("Scar cells:", int(scar_mask.sum()), "out of", full.n_cells)

    output_path = write_vtk_ascii(full, options.output_path)

    postprocess_vtk_output(
        output_path,
        convert_fields=options.convert_fields,
        remove_blanks=options.remove_blank_lines,
    )
    if options.inspect:
        inspect_vtk(output_path, "Inspecting Output Fields")
    print(f"Scar-tagged mesh written to {output_path}")
    return StepResult(output_mesh=output_path)


def run_from_context(context: StepContext) -> StepResult:
    opts = context.options
    options = ScarOptions(
        full_mesh_path=opts.get("full_mesh_path") or (context.current_mesh or ""),
        selection_path=opts["selection_path"],
        output_path=opts["output_path"],
        id_array=str(opts.get("id_array", "GlobalCellIds")),
        scar_name=str(opts.get("scar_name", "Scar")),
        scar_value=float(opts.get("scar_value", 1.0)),
        diffusivity_name=str(opts.get("diffusivity_name", "Diffusivity")),
        diffusivity_scale=float(opts.get("diffusivity_scale", 0.0)),
        inspect=bool(opts.get("inspect", False)),
        convert_fields=bool(opts.get("convert_fields", True)),
        remove_blank_lines=bool(opts.get("remove_blank_lines", True)),
    )
    return run_scar(options)
