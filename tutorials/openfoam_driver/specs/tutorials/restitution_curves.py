from __future__ import annotations

import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path

from ...core.defaults import restitution_curves as defaults
from ...postprocessing.driver import PostprocessTask, run_postprocess_tasks
from ..common import (
    resolve_run_script_path,
    resolve_spec_paths,
    set_ionic_model,
    set_n_stim1,
    set_n_stim2,
    set_s1_period,
    set_s2_period,
    set_stimulus_amplitude,
    set_tissue,
    set_write_after_time,
    set_end_time,
)
from ...core.runtime.models import CaseConfig, TutorialSpec


def _build_cases(
    ionic_models: Sequence[str],
    ionic_model_tissue_map: Mapping[str, Sequence[str]],
    s2_intervals_ms: Sequence[int],
) -> list[CaseConfig]:
    cases: list[CaseConfig] = []
    for ionic_model in ionic_models:
        tissues = ionic_model_tissue_map.get(ionic_model)
        if not tissues:
            raise KeyError(f"Missing tissue list for ionic model '{ionic_model}'")
        for tissue in tissues:
            for s2 in s2_intervals_ms:
                case_id = f"{ionic_model}_{tissue}_S2_{s2}"
                cases.append(
                    CaseConfig(
                        case_id=case_id,
                        params={
                            "ionicModel": ionic_model,
                            "tissue": tissue,
                            "s2Interval": s2,
                        },
                    )
                )
    return cases


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    stimulus_map: Mapping[str, float],
    s1_interval_ms: int,
    n_s1: int,
    n_s2: int,
    write_after_time_s: float,
    end_time_buffer_s: float,
    properties_relpath: Path = defaults.CARDIAC_PROPERTIES_RELPATH,
    stimulus_relpath: Path = defaults.STIMULUS_RELPATH,
    control_dict_relpath: Path = defaults.CONTROL_DICT_RELPATH,
) -> None:
    ionic_model = case.params["ionicModel"]
    tissue = case.params["tissue"]
    s2_interval_ms = case.params["s2Interval"]

    if ionic_model not in stimulus_map:
        raise KeyError(f"Missing stimulus amplitude for ionic model '{ionic_model}'")

    properties_file = case_root / properties_relpath
    stimulus_file = case_root / stimulus_relpath
    control_dict_file = case_root / control_dict_relpath

    # cardiacProperties
    set_tissue(properties_file, tissue)
    set_ionic_model(properties_file, ionic_model)

    # stimulusProtocol
    set_stimulus_amplitude(stimulus_file, stimulus_map[ionic_model])
    set_s1_period(stimulus_file, s1_interval_ms)
    set_n_stim1(stimulus_file, n_s1)
    set_s2_period(stimulus_file, s2_interval_ms)
    set_n_stim2(stimulus_file, n_s2)
    set_write_after_time(stimulus_file, write_after_time_s)

    # controlDict: endTime = (S1*n_S1 + S2*n_S2) / 1000 + buffer
    end_time = (s1_interval_ms * n_s1 + s2_interval_ms * n_s2) / 1000.0 + end_time_buffer_s
    set_end_time(control_dict_file, end_time)


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    tutorials_root: Path | None = None,
    run_script_relpath: Path = defaults.RUN_SCRIPT_RELPATH,
    output_dir: Path,
) -> None:
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )
    subprocess.run(
        ["bash", "-l", str(run_script), "--case-dir", str(case_root)],
        check=True,
    )

    # Per-case collection: generate video, then rename & move the .txt output
    # before the next simulation overwrites it.
    ionic_model = case.params["ionicModel"]
    tissue = case.params["tissue"]
    s2_interval = case.params["s2Interval"]
    output_dir.mkdir(parents=True, exist_ok=True)

    for source in sorted(case_root.glob("*.txt")):
        stem = f"{ionic_model}_{tissue}_S2_{s2_interval}"

        # Animate the voltage trace before moving the file
        _animate_voltage_trace(
            data_file=source,
            output_path=output_dir / f"{stem}.mp4",
            title=f"{ionic_model} · {tissue} · S2={s2_interval} ms",
        )

        dest = output_dir / f"{stem}.txt"
        if dest.exists():
            dest.unlink()
        shutil.move(str(source), str(dest))
        print(f"Moved output: {source.name} -> {dest}")


def _animate_voltage_trace(
    data_file: Path,
    output_path: Path,
    title: str = "",
) -> None:
    """Load animate_trace.py from the tutorial setup dir and call it."""
    # The module lives in the tutorial's setupRestitutionCurves_s1s2Protocol/
    # directory, next to this spec's sibling scripts. We discover it relative
    # to the data file's grandparent (case_root).
    animate_module = data_file.parent.parent / "setupRestitutionCurves_s1s2Protocol" / "animate_trace.py"
    if not animate_module.exists():
        print(f"Warning: animate_trace.py not found at {animate_module}; skipping video.")
        return

    import importlib.util
    mod_spec = importlib.util.spec_from_file_location("animate_trace", animate_module)
    if mod_spec is None or mod_spec.loader is None:
        print(f"Warning: could not load animate_trace module; skipping video.")
        return
    module = importlib.util.module_from_spec(mod_spec)
    mod_spec.loader.exec_module(module)
    module.create_voltage_animation(data_file, output_path, title=title)


def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    postprocess_script_relpath: Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    ionic_models: Sequence[str],
    ionic_model_tissue_map: Mapping[str, Sequence[str]],
    show_plots: bool = False,
    strict_artifacts: bool = False,
) -> None:
    run_postprocess_tasks(
        setup_root=setup_root,
        output_dir=output_dir,
        tutorial_name=defaults.TUTORIAL_NAME,
        strict_artifacts=strict_artifacts,
        tasks=[
            PostprocessTask(
                module_relpath=postprocess_script_relpath,
                function_name=postprocess_function_name,
                kwargs={
                    "ionic_models": list(ionic_models),
                    "tissue_map": {m: list(t) for m, t in ionic_model_tissue_map.items()},
                    "show_plots": show_plots,
                },
            )
        ],
    )


def make_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str = defaults.TUTORIAL_NAME,
    setup_dir_name: str | None = None,
    output_dir_name: str | None = None,
    ionic_models: Sequence[str] = defaults.IONIC_MODELS,
    ionic_model_tissue_map: Mapping[str, Sequence[str]] = defaults.IONIC_MODEL_TISSUE_MAP,
    stimulus_map: Mapping[str, float] = defaults.STIMULUS_MAP,
    s1_interval_ms: int = defaults.S1_INTERVAL_MS,
    n_s1: int = defaults.N_S1,
    n_s2: int = defaults.N_S2,
    s2_intervals_ms: Sequence[int] = defaults.S2_INTERVALS_MS,
    end_time_buffer_s: float = defaults.END_TIME_BUFFER_S,
    properties_relpath: str | Path = defaults.CARDIAC_PROPERTIES_RELPATH,
    stimulus_relpath: str | Path = defaults.STIMULUS_RELPATH,
    control_dict_relpath: str | Path = defaults.CONTROL_DICT_RELPATH,
    run_script_relpath: str | Path = defaults.RUN_SCRIPT_RELPATH,
    output_glob: str = defaults.OUTPUT_GLOB,
    postprocess_script_relpath: str | Path = defaults.POSTPROCESS_SCRIPT_RELPATH,
    postprocess_function_name: str = defaults.POSTPROCESS_FUNCTION_NAME,
    show_plots: bool = False,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    ionic_models_list = [str(m) for m in ionic_models]
    if not ionic_models_list:
        raise ValueError("ionic_models cannot be empty")

    for ionic_model in ionic_models_list:
        if not ionic_model_tissue_map.get(ionic_model):
            raise KeyError(f"Missing tissue mapping for ionic model '{ionic_model}'")
        if ionic_model not in stimulus_map:
            raise KeyError(f"Missing stimulus amplitude for ionic model '{ionic_model}'")

    if not s2_intervals_ms:
        raise ValueError("s2_intervals_ms cannot be empty")

    properties_path = Path(properties_relpath)
    stimulus_path = Path(stimulus_relpath)
    control_dict_path = Path(control_dict_relpath)
    run_script_path = Path(run_script_relpath)
    postprocess_script_path = Path(postprocess_script_relpath)

    # write_after_time: start writing 2 s before the end of the S1 phase
    write_after_time_s = (s1_interval_ms * n_s1) / 1000.0 - 2.0

    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=defaults.OUTPUT_DIR_NAME,
    )

    return TutorialSpec(
        name=case_dir_name,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=partial(
            _build_cases,
            ionic_models=ionic_models_list,
            ionic_model_tissue_map=ionic_model_tissue_map,
            s2_intervals_ms=list(s2_intervals_ms),
        ),
        apply_case=partial(
            _apply_case,
            stimulus_map=stimulus_map,
            s1_interval_ms=s1_interval_ms,
            n_s1=n_s1,
            n_s2=n_s2,
            write_after_time_s=write_after_time_s,
            end_time_buffer_s=end_time_buffer_s,
            properties_relpath=properties_path,
            stimulus_relpath=stimulus_path,
            control_dict_relpath=control_dict_path,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=run_script_path,
            output_dir=output_dir,
        ),
        collect_outputs=None,  # per-case collection handled inside _run_case
        postprocess=partial(
            _postprocess,
            postprocess_script_relpath=postprocess_script_path,
            postprocess_function_name=postprocess_function_name,
            ionic_models=ionic_models_list,
            ionic_model_tissue_map=ionic_model_tissue_map,
            show_plots=show_plots,
            strict_artifacts=postprocess_strict_artifacts,
        ),
        metadata={
            "python": sys.executable,
            "notes": "S1–S2 restitution protocol sweep on ionic model, tissue, and S2 interval.",
            "ionic_models": ionic_models_list,
            "s1_interval_ms": s1_interval_ms,
            "n_s1": n_s1,
            "n_s2": n_s2,
            "s2_intervals_ms": list(s2_intervals_ms),
            "properties_relpath": str(properties_path),
            "stimulus_relpath": str(stimulus_path),
            "control_dict_relpath": str(control_dict_path),
            "run_script_relpath": str(run_script_path),
            "output_glob": output_glob,
            "postprocess_script_relpath": str(postprocess_script_path),
            "postprocess_function_name": postprocess_function_name,
            "show_plots": show_plots,
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
