from __future__ import annotations

import subprocess
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path
from typing import Any

from ...core.defaults.shared import OUTPUT_DIR_NAME, RUN_CASE_SCRIPT_RELPATH
from ...core.runtime.models import CaseConfig, TutorialSpec
from ...postprocessing.driver import PostprocessTask, run_postprocess_tasks
from ..common import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
    collect_outputs_by_pattern,
    resolve_run_script_path,
    resolve_spec_paths,
)


def _normalize_case_specs(
    *,
    cases: Sequence[Mapping[str, Any]] | None,
    electro_property_overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
    physics_property_overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
    dimension: str | None,
    parallel: bool,
    touch_case_foam: bool,
    openfoam_bashrc: str | Path | None,
) -> list[CaseConfig]:
    if cases is None:
        payload = {
            "electro_property_overrides": electro_property_overrides,
            "physics_property_overrides": physics_property_overrides,
            "dimension": dimension,
            "parallel": parallel,
            "touch_case_foam": touch_case_foam,
            "openfoam_bashrc": str(openfoam_bashrc) if openfoam_bashrc is not None else None,
        }
        return [CaseConfig(case_id="default", params=payload)]

    normalized: list[CaseConfig] = []
    for index, item in enumerate(cases, start=1):
        case_id = str(item.get("case_id", f"case{index:03d}"))
        normalized.append(
            CaseConfig(
                case_id=case_id,
                params={
                    "electro_property_overrides": item.get(
                        "electro_property_overrides",
                        electro_property_overrides,
                    ),
                    "physics_property_overrides": item.get(
                        "physics_property_overrides",
                        physics_property_overrides,
                    ),
                    "dimension": item.get("dimension", dimension),
                    "parallel": bool(item.get("parallel", parallel)),
                    "touch_case_foam": bool(item.get("touch_case_foam", touch_case_foam)),
                    "openfoam_bashrc": (
                        str(item["openfoam_bashrc"])
                        if item.get("openfoam_bashrc") is not None
                        else (str(openfoam_bashrc) if openfoam_bashrc is not None else None)
                    ),
                },
            )
        )
    return normalized


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    electro_properties_relpath: Path,
    physics_properties_relpath: Path,
) -> None:
    apply_electro_property_overrides(
        case_root / electro_properties_relpath,
        case.params.get("electro_property_overrides"),
    )
    apply_physics_property_overrides(
        case_root / physics_properties_relpath,
        case.params.get("physics_property_overrides"),
    )


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    tutorials_root: Path | None,
    run_script_relpath: Path,
) -> None:
    del setup_root
    run_script = resolve_run_script_path(
        tutorials_root=tutorials_root,
        run_script_relpath=run_script_relpath,
    )

    command = [
        "bash",
        "-l",
        str(run_script),
        "--case-dir",
        str(case_root),
    ]

    dimension = case.params.get("dimension")
    if dimension:
        command.extend(["--dimension", str(dimension)])
    if case.params.get("parallel"):
        command.append("--parallel")
    if case.params.get("touch_case_foam"):
        command.append("--touch-case-foam")

    openfoam_bashrc = case.params.get("openfoam_bashrc")
    if openfoam_bashrc:
        command.extend(["--openfoam-bashrc", str(openfoam_bashrc)])

    subprocess.run(command, check=True)


def _collect_outputs(case_root: Path, output_dir: Path, *, patterns: Sequence[str]) -> None:
    for pattern in patterns:
        collect_outputs_by_pattern(case_root, output_dir, pattern=pattern)


def _postprocess(
    setup_root: Path,
    output_dir: Path,
    *,
    tutorial_name: str,
    postprocess_tasks: Sequence[PostprocessTask],
    strict_artifacts: bool,
) -> None:
    run_postprocess_tasks(
        setup_root=setup_root,
        output_dir=output_dir,
        tutorial_name=tutorial_name,
        tasks=list(postprocess_tasks),
        strict_artifacts=strict_artifacts,
    )


def _normalize_postprocess_tasks(
    tasks: Sequence[Mapping[str, Any]] | None,
) -> list[PostprocessTask]:
    normalized: list[PostprocessTask] = []
    for item in tasks or ():
        if "module_relpath" not in item:
            raise KeyError("Generic postprocess tasks require 'module_relpath'")
        normalized.append(
            PostprocessTask(
                module_relpath=Path(str(item["module_relpath"])),
                function_name=str(item.get("function_name", "run_postprocessing")),
                kwargs=dict(item.get("kwargs", {})),
            )
        )
    return normalized


def make_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str,
    setup_dir_name: str | None = None,
    output_dir_name: str | None = None,
    electro_properties_relpath: str | Path = "constant/electroProperties",
    physics_properties_relpath: str | Path = "constant/physicsProperties",
    electro_property_overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None = None,
    physics_property_overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None = None,
    cases: Sequence[Mapping[str, Any]] | None = None,
    dimension: str | None = None,
    parallel: bool = False,
    touch_case_foam: bool = False,
    openfoam_bashrc: str | Path | None = None,
    collect_patterns: Sequence[str] = (),
    postprocess_tasks: Sequence[Mapping[str, Any]] | None = None,
    run_script_relpath: str | Path = RUN_CASE_SCRIPT_RELPATH,
    postprocess_strict_artifacts: bool = False,
) -> TutorialSpec:
    if not str(case_dir_name).strip():
        raise ValueError("case_dir_name cannot be empty")

    electro_properties_path = Path(electro_properties_relpath)
    physics_properties_path = Path(physics_properties_relpath)
    run_script_path = Path(run_script_relpath)
    normalized_postprocess_tasks = _normalize_postprocess_tasks(postprocess_tasks)

    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=OUTPUT_DIR_NAME,
    )

    normalized_cases = _normalize_case_specs(
        cases=cases,
        electro_property_overrides=electro_property_overrides,
        physics_property_overrides=physics_property_overrides,
        dimension=dimension,
        parallel=parallel,
        touch_case_foam=touch_case_foam,
        openfoam_bashrc=openfoam_bashrc,
    )

    return TutorialSpec(
        name=case_dir_name,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=lambda: list(normalized_cases),
        apply_case=partial(
            _apply_case,
            electro_properties_relpath=electro_properties_path,
            physics_properties_relpath=physics_properties_path,
        ),
        run_case=partial(
            _run_case,
            tutorials_root=tutorials_root,
            run_script_relpath=run_script_path,
        ),
        collect_outputs=(
            partial(_collect_outputs, patterns=tuple(str(item) for item in collect_patterns))
            if collect_patterns
            else None
        ),
        postprocess=(
            partial(
                _postprocess,
                tutorial_name=case_dir_name,
                postprocess_tasks=tuple(normalized_postprocess_tasks),
                strict_artifacts=postprocess_strict_artifacts,
            )
            if normalized_postprocess_tasks
            else None
        ),
        metadata={
            "notes": "Generic case runner for arbitrary tutorial folders.",
            "electro_properties_relpath": str(electro_properties_path),
            "physics_properties_relpath": str(physics_properties_path),
            "run_script_relpath": str(run_script_path),
            "collect_patterns": list(collect_patterns),
            "postprocess_task_count": len(normalized_postprocess_tasks),
            "case_count": len(normalized_cases),
            "has_default_electro_property_overrides": bool(electro_property_overrides),
            "has_default_physics_property_overrides": bool(physics_property_overrides),
            "postprocess_strict_artifacts": postprocess_strict_artifacts,
        },
    )
