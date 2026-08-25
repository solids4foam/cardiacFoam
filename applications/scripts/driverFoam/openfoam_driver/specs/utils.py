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


def set_delta_t(control_dict_path: Path, delta_t_seconds: float) -> None:
    update_foam_entry(control_dict_path, "deltaT", delta_t_seconds)


def set_end_time(control_dict_path: Path, t_s: float) -> None:
    update_foam_entry(control_dict_path, "endTime", t_s)



def replace_block_mesh_resolutions(
    block_mesh_dict_path: Path,
    cell_counts_str: str,
    *,
    expected_blocks: int = 1,
) -> None:
    """Rewrite lines starting with ``hex (`` in an existing `block_mesh_dict_path`.
    Replaces the cell counts portion of the hex definition with `cell_counts_str`.
    Validates that exactly `expected_blocks` were replaced.
    """
    if not block_mesh_dict_path.exists():
        raise FileNotFoundError(f"Missing mesh dictionary: {block_mesh_dict_path}")

    lines = block_mesh_dict_path.read_text().splitlines(keepends=True)
    replaced_count = 0

    with block_mesh_dict_path.open("w") as handle:
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("hex (") and not stripped.startswith("//"):
                prefix, _, suffix = line.partition(") (")
                if not suffix:
                    handle.write(line)
                    continue
                _, _, trailing = suffix.partition(") simpleGrading")
                handle.write(f"{prefix}) ({cell_counts_str}) simpleGrading{trailing}")
                replaced_count += 1
            else:
                handle.write(line)

    if replaced_count != expected_blocks:
        raise KeyError(
            f"Expected to update {expected_blocks} hex blocks in {block_mesh_dict_path}, "
            f"but found {replaced_count}."
        )


def archive_case_logs(case_root: Path, case_id: str) -> Path | None:
    log_files = sorted(path for path in case_root.glob("log.*") if path.is_file())
    if not log_files:
        return None

    destination_root = case_root / "logs" / case_id
    if destination_root.exists():
        shutil.rmtree(destination_root)
    destination_root.mkdir(parents=True, exist_ok=True)

    for source in log_files:
        shutil.copy2(source, destination_root / source.name)

    print(f"Archived {len(log_files)} log file(s) for {case_id}: {destination_root}")
    return destination_root


def stage_post_processing_outputs(
    case_root: Path,
    destination_dir: Path,
    file_mapping: dict[str, str],
    *,
    missing_ok: bool = False,
) -> list[Path]:
    """
    Finds files in `case_root/postProcessing/` (or `processor0/postProcessing/` fallback)
    and copies or moves them to `destination_dir` using the renamed target names.
    file_mapping is a dict of {source_filename: destination_filename}.
    """
    staged_outputs: list[Path] = []
    destination_dir.mkdir(parents=True, exist_ok=True)

    for source_name, dest_name in file_mapping.items():
        destination = destination_dir / dest_name
        candidates = (
            case_root / "postProcessing" / source_name,
            case_root / "processor0" / "postProcessing" / source_name,
        )
        found = False
        for candidate in candidates:
            if not candidate.exists():
                continue
            if candidate.parent == destination_dir:
                shutil.move(str(candidate), str(destination))
            else:
                shutil.copy2(candidate, destination)
            print(f"Archived output: {candidate} -> {destination}")
            staged_outputs.append(destination)
            found = True
            break

        if not found and not missing_ok:
            checked = ", ".join(str(path) for path in candidates)
            raise FileNotFoundError(
                f"Output '{source_name}' not found after run. Checked: {checked}"
            )

    return staged_outputs

