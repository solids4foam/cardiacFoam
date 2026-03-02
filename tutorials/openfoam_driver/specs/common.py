from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

from ..core.runtime.mutators import update_foam_entry


def tutorials_root_default() -> Path:
    return Path(__file__).resolve().parents[2]


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
    resolved_tutorials_root = tutorials_root or tutorials_root_default()
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

    resolved_tutorials_root = tutorials_root or tutorials_root_default()
    driver_level_path = resolved_tutorials_root / run_script_relpath
    if not driver_level_path.exists():
        raise FileNotFoundError(
            f"Run script not found: {driver_level_path}. "
            "Use an absolute path or a path relative to tutorials root."
        )
    return driver_level_path


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


def set_cardiac_dimension(cardiac_properties_path: Path, dimension: str) -> None:
    update_foam_entry(cardiac_properties_path, "dimension", f'"{dimension}"')


def set_tissue(cardiac_properties_path: Path, tissue: str) -> None:
    update_foam_entry(cardiac_properties_path, "tissue", tissue)


def set_ionic_model(cardiac_properties_path: Path, ionic_model: str) -> None:
    update_foam_entry(cardiac_properties_path, "ionicModel", ionic_model)


def set_stimulus_amplitude(stimulus_file_path: Path, amplitude: float) -> None:
    update_foam_entry(stimulus_file_path, "stim_amplitude", amplitude)


def solver_to_explicit_flag(solver: str) -> str:
    return "yes" if solver == "explicit" else "no"


def set_solve_explicit(time_integration_properties_path: Path, solver: str) -> None:
    update_foam_entry(
        time_integration_properties_path,
        "solveExplicit",
        solver_to_explicit_flag(solver),
    )


def set_s1_period(stimulus_path: Path, s1_interval_ms: int | float) -> None:
    update_foam_entry(stimulus_path, "stim_period_S1", s1_interval_ms)


def set_s2_period(stimulus_path: Path, s2_interval_ms: int | float) -> None:
    update_foam_entry(stimulus_path, "stim_period_S2", s2_interval_ms)


def set_n_stim1(stimulus_path: Path, n: int) -> None:
    update_foam_entry(stimulus_path, "nstim1", n)


def set_n_stim2(stimulus_path: Path, n: int) -> None:
    update_foam_entry(stimulus_path, "nstim2", n)


def set_write_after_time(stimulus_path: Path, t_s: float) -> None:
    update_foam_entry(stimulus_path, "writeAfterTime", t_s)


def set_end_time(control_dict_path: Path, t_s: float) -> None:
    update_foam_entry(control_dict_path, "endTime", t_s)
