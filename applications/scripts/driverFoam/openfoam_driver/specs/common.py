from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

from ..core.runtime.mutators import update_foam_entry


def repo_root_default() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tutorials").exists() and (parent / "src").exists():
            return parent
    return current.parents[2]


def tutorials_root_default() -> Path:
    repo_root = repo_root_default()
    tutorials_root = repo_root / "tutorials"
    if tutorials_root.exists():
        return tutorials_root
    return repo_root


def default_setup_dir_name(case_dir_name: str) -> str:
    normalized_case_dir = case_dir_name.strip()
    if not normalized_case_dir:
        raise ValueError("case_dir_name cannot be empty")
    return f"setup{normalized_case_dir[:1].upper()}{normalized_case_dir[1:]}"


def resolve_spec_paths(
    *,
    tutorials_root: Path | None,
    case_dir_name: str,
    setup_dir_name: str | None = None,
    output_dir_name: str | Path | None = None,
    default_output_dir_name: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    resolved_tutorials_root = (
        Path(tutorials_root) if tutorials_root is not None else tutorials_root_default()
    )
    resolved_case_dir = case_dir_name.strip()
    if not resolved_case_dir:
        raise ValueError("case_dir_name cannot be empty")
    resolved_setup_dir = setup_dir_name or default_setup_dir_name(resolved_case_dir)
    resolved_output_dir = output_dir_name or default_output_dir_name
    if resolved_output_dir is None:
        raise ValueError("output_dir_name or default_output_dir_name must be provided")

    case_root = resolved_tutorials_root / resolved_case_dir
    setup_root = case_root / resolved_setup_dir
    output_dir = case_root / Path(resolved_output_dir)
    return case_root, setup_root, output_dir


def resolve_run_script_path(
    *,
    tutorials_root: Path | None,
    run_script_relpath: Path,
) -> Path:
    if run_script_relpath.is_absolute():
        return run_script_relpath

    candidate_roots: list[Path] = []
    if tutorials_root is not None:
        candidate_roots.append(Path(tutorials_root))
    candidate_roots.append(repo_root_default())
    candidate_roots.append(tutorials_root_default())

    checked_paths: list[Path] = []
    seen: set[Path] = set()
    for root in candidate_roots:
        resolved_root = root.resolve()
        if resolved_root in seen:
            continue
        seen.add(resolved_root)
        candidate = resolved_root / run_script_relpath
        checked_paths.append(candidate)
        if candidate.exists():
            return candidate

    checked_str = ", ".join(str(path) for path in checked_paths)
    raise FileNotFoundError(
        f"Run script not found for '{run_script_relpath}'. Checked: {checked_str}. "
        "Use an absolute path or a path relative to repository root."
    )


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


def set_cardiac_dimension(
    electro_properties_path: Path,
    dimension: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(
        electro_properties_path,
        "dimension",
        f'"{dimension}"',
        scope=scope,
    )


def set_tissue(
    electro_properties_path: Path,
    tissue: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(electro_properties_path, "tissue", tissue, scope=scope)


def set_ionic_model(
    electro_properties_path: Path,
    ionic_model: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(
        electro_properties_path,
        "ionicModel",
        ionic_model,
        scope=scope,
    )


def set_stimulus_amplitude(
    electro_properties_path: Path,
    amplitude: float,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(
        electro_properties_path,
        "stim_amplitude",
        amplitude,
        scope=scope,
    )


def set_solution_algorithm(
    electro_properties_path: Path,
    solver: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    algorithm = "explicit" if solver == "explicit" else "implicit"
    update_foam_entry(
        electro_properties_path,
        "solutionAlgorithm",
        algorithm,
        scope=scope,
    )


def set_s1_period(
    electro_properties_path: Path,
    s1_interval_ms: int | float,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(
        electro_properties_path,
        "stim_period_S1",
        s1_interval_ms,
        scope=scope,
    )


def set_s2_period(
    electro_properties_path: Path,
    s2_interval_ms: int | float,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(
        electro_properties_path,
        "stim_period_S2",
        s2_interval_ms,
        scope=scope,
    )


def set_n_stim1(
    electro_properties_path: Path,
    n: int,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(electro_properties_path, "nstim1", n, scope=scope)


def set_n_stim2(
    electro_properties_path: Path,
    n: int,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(electro_properties_path, "nstim2", n, scope=scope)


def set_write_after_time(
    electro_properties_path: Path,
    t_s: float,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    update_foam_entry(electro_properties_path, "writeAfterTime", t_s, scope=scope)


def set_end_time(control_dict_path: Path, t_s: float) -> None:
    update_foam_entry(control_dict_path, "endTime", t_s)
