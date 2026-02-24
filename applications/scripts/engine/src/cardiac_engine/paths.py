"""Canonical repository paths for Phase 1 architecture."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RepoPaths:
    scripts_root: Path
    project_root: Path
    files_organize: Path

    diffusivity_config: Path
    scar_config: Path
    purkinje_slab_config: Path
    purkinje_fractal_config: Path
    purkinje_fractal_runner: Path
    convert_config: Path
    convert_runner: Path

    default_input_mesh: Path
    default_fractal_mesh: Path
    default_outputs_root: Path


def resolve_paths(*, project_root: Path) -> RepoPaths:
    scripts_root = project_root.parent
    files_organize = scripts_root / "filesOrganize"
    return RepoPaths(
        scripts_root=scripts_root,
        project_root=project_root,
        files_organize=files_organize,
        diffusivity_config=scripts_root / "diffusivity" / "Diffusivity_config.py",
        scar_config=scripts_root / "scar" / "scar_config.py",
        purkinje_slab_config=scripts_root / "purkinje" / "slab" / "purkinjeSlab_config.py",
        purkinje_fractal_config=scripts_root / "purkinje" / "fractal_3d" / "purkinjeFractalTree_config.py",
        purkinje_fractal_runner=scripts_root / "purkinje" / "fractal_3d" / "purkinje3D_fractal.py",
        convert_config=project_root / "configs" / "convert_config.py",
        convert_runner=project_root / "src" / "cardiac_preproc" / "fileConversion" / "ASCIIlegacyToVtkUnstructured.py",
        default_input_mesh=files_organize / "cardiac_preproc" / "inputs" / "meshes" / "ASCIIlegacy.vtk",
        default_fractal_mesh=files_organize / "cardiac_preproc" / "inputs" / "meshes" / "biv_ellipsoid.msh",
        default_outputs_root=files_organize / "cardiac_preproc" / "outputs",
    )
