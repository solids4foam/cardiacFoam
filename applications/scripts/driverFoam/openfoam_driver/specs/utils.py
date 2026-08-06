import importlib.util
import shutil
from collections.abc import Mapping
from pathlib import Path

from ..core.runtime.mutators import update_foam_entry


def load_python_module(module_path: Path, *, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collect_outputs_by_pattern(case_root: Path, output_dir: Path, *, pattern: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    matching_files = sorted(case_root.glob(pattern))
    for source in matching_files:
        destination = output_dir / source.name
        if destination.exists():
            destination.unlink()
        shutil.move(str(source), str(destination))
        print(f"Moved output: {source.name} -> {destination}")


def set_delta_t(control_dict_path: Path, delta_t_seconds: float) -> None:
    update_foam_entry(control_dict_path, "deltaT", delta_t_seconds)


def set_end_time(control_dict_path: Path, t_s: float) -> None:
    update_foam_entry(control_dict_path, "endTime", t_s)


def replace_single_block_mesh_resolution(
    block_mesh_dict_path: Path,
    cells: int,
    dimension: str,
    *,
    resolution_by_dimension: Mapping[str, str],
) -> None:
    """Rewrite the single ``hex (0 1 2 3 4 5 6 7) (...)`` line in an existing
    `block_mesh_dict_path` via a cells+dimension template lookup.
    """
    if not block_mesh_dict_path.exists():
        raise FileNotFoundError(f"Missing mesh dictionary: {block_mesh_dict_path}")

    try:
        cell_counts = resolution_by_dimension[dimension].format(cells=cells)
    except KeyError as exc:
        raise ValueError(f"Unsupported dimension: {dimension}") from exc

    replacement = f"hex (0 1 2 3 4 5 6 7) ({cell_counts}) simpleGrading (1 1 1)\n"

    lines = block_mesh_dict_path.read_text().splitlines(keepends=True)
    replaced = False
    with block_mesh_dict_path.open("w") as handle:
        for line in lines:
            stripped = line.strip()
            if (
                not replaced
                and stripped.startswith("hex (0 1 2 3 4 5 6 7)")
                and not stripped.startswith("//")
            ):
                handle.write(replacement)
                replaced = True
            else:
                handle.write(line)

    if not replaced:
        raise KeyError(f"Target hex line not found in {block_mesh_dict_path}")
