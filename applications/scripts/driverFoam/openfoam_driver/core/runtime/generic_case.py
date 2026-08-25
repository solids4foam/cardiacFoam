"""Solver-neutral generic case spec implementation.

This module owns the default case-folder execution contract used by the
registry for arbitrary OpenFOAM folders. Plugins may still call the same
factory for richer flows such as build-and-launch, but the implementation
itself lives in core.
"""

from __future__ import annotations

import subprocess
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path
from typing import Any

from openfoam_driver.core.specs.common import (
    resolve_run_script_path,
    resolve_spec_paths,
)

from .models import CaseConfig, TutorialSpec

OUTPUT_DIR_NAME = "postProcessing"
RUN_CASE_SCRIPT_RELPATH = Path(
    "applications/scripts/driverFoam/openfoam_driver/scripts/run_case.sh"
)


DictFileOverrides = Mapping[str, "Mapping[str, Any] | Sequence[Mapping[str, Any]]"]


def _merged_dict_file_overrides(
    item: Mapping[str, Any],
    default_overrides: Mapping[str, Any],
) -> dict[str, Any]:
    """Per-case ``dict_file_overrides``, falling back to the spec-level default.

    A case entry may name only some of the dictionary files; the remaining ones
    keep whatever the spec-level default declared for them.
    """
    from openfoam_driver.core.compatibility import (
        legacy_generic_case_alias_names,
        legacy_generic_case_dict_file_aliases,
    )

    merged = dict(default_overrides)
    merged.update(item.get("dict_file_overrides") or {})

    deprecated = {
        key: value
        for key, value in item.items()
        if key in legacy_generic_case_alias_names()
    }
    if deprecated:
        _relpaths, alias_overrides, _unknown = legacy_generic_case_dict_file_aliases(
            deprecated,
        )
        merged.update(alias_overrides)
    return {key: value for key, value in merged.items() if value}


def _normalize_case_specs(
    *,
    cases: Sequence[Mapping[str, Any]] | None,
    dict_file_overrides: Mapping[str, Any],
    dimension: str | None,
    parallel: bool,
    touch_case_foam: bool,
    openfoam_bashrc: str | Path | None,
    solver_command: str | Sequence[str] | None,
    pre_solve_commands: Sequence[str | Sequence[str]],
) -> list[CaseConfig]:
    if cases is None:
        payload = {
            "dict_file_overrides": dict(dict_file_overrides),
            "dimension": dimension,
            "parallel": parallel,
            "touch_case_foam": touch_case_foam,
            "openfoam_bashrc": str(openfoam_bashrc) if openfoam_bashrc is not None else None,
            "solver_command": solver_command,
            "pre_solve_commands": list(pre_solve_commands),
        }
        return [CaseConfig(case_id="default", params=payload)]

    normalized: list[CaseConfig] = []
    for index, item in enumerate(cases, start=1):
        case_id = str(item.get("case_id", f"case{index:03d}"))
        normalized.append(
            CaseConfig(
                case_id=case_id,
                params={
                    "dict_file_overrides": _merged_dict_file_overrides(
                        item, dict_file_overrides,
                    ),
                    "dimension": item.get("dimension", dimension),
                    "parallel": bool(item.get("parallel", parallel)),
                    "touch_case_foam": bool(item.get("touch_case_foam", touch_case_foam)),
                    "openfoam_bashrc": (
                        str(item["openfoam_bashrc"])
                        if item.get("openfoam_bashrc") is not None
                        else (str(openfoam_bashrc) if openfoam_bashrc is not None else None)
                    ),
                    "solver_command": item.get("solver_command", solver_command),
                    "pre_solve_commands": list(item.get("pre_solve_commands", pre_solve_commands)),
                },
            )
        )
    return normalized


def _apply_case(
    case_root: Path,
    case: CaseConfig,
    *,
    dict_file_relpaths: Mapping[str, Path],
    mutation_callback,
) -> None:
    mutation_callback(
        case_root,
        case,
        dict_file_relpaths=dict(dict_file_relpaths),
        dict_file_overrides=dict(case.params.get("dict_file_overrides") or {}),
    )


def _split_command(command: str | Sequence[str]) -> list[str]:
    return command.split() if isinstance(command, str) else list(command)


def _workflow_dag_for(
    *,
    solver_command: str | Sequence[str] | None,
    pre_solve_commands: Sequence[str | Sequence[str]],
) -> dict[str, Any]:
    if solver_command is None:
        return {"steps": [{"id": "run", "command": "Allrun", "depends_on": []}]}

    steps: list[dict[str, Any]] = []
    depends_on: list[str] = []
    for index, raw_cmd in enumerate(pre_solve_commands):
        cmd = _split_command(raw_cmd)
        step_id = f"pre_{index}"
        steps.append(
            {
                "id": step_id,
                "command": cmd[0],
                "args": cmd[1:],
                "depends_on": depends_on,
            }
        )
        depends_on = [step_id]

    solve_cmd = _split_command(solver_command)
    steps.append(
        {
            "id": "solve",
            "command": solve_cmd[0],
            "args": solve_cmd[1:],
            "depends_on": depends_on,
        }
    )
    return {"steps": steps}


def make_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str,
    setup_dir_name: str | None = None,
    output_dir_name: str | None = None,
    dict_file_relpaths: Mapping[str, str | Path] | None = None,
    dict_file_overrides: DictFileOverrides | None = None,
    cases: Sequence[Mapping[str, Any]] | None = None,
    dimension: str | None = None,
    parallel: bool = False,
    touch_case_foam: bool = False,
    openfoam_bashrc: str | Path | None = None,
    collect_patterns: Sequence[str] = (),
    run_script_relpath: str | Path = RUN_CASE_SCRIPT_RELPATH,
    solver_command: str | Sequence[str] | None = None,
    pre_solve_commands: Sequence[str | Sequence[str]] | None = None,
    _apply_case_mutation=None,
    **deprecated_kwargs: Any,
) -> TutorialSpec:
    if not str(case_dir_name).strip():
        raise ValueError("case_dir_name cannot be empty")

    from openfoam_driver.core.compatibility import (
        legacy_generic_case_dict_file_aliases,
        legacy_generic_case_dict_file_relpaths,
    )

    # ``dict_file_relpaths``/``dict_file_overrides`` are the supported spelling.
    # The historical solver-specific keyword names survive only as deprecated
    # aliases, translated by the named compatibility seam so that no plugin
    # vocabulary reaches core.
    alias_relpaths: dict[str, Any] = {}
    alias_overrides: dict[str, Any] = {}
    unknown_kwargs: list[str] = []
    if deprecated_kwargs:
        alias_relpaths, alias_overrides, unknown_kwargs = (
            legacy_generic_case_dict_file_aliases(deprecated_kwargs)
        )
    if unknown_kwargs:
        raise TypeError(
            "make_spec() got an unexpected keyword argument "
            + ", ".join(repr(key) for key in sorted(unknown_kwargs))
        )

    resolved_relpaths_raw: dict[str, Any] = dict(
        legacy_generic_case_dict_file_relpaths()
        if dict_file_relpaths is None
        else dict_file_relpaths
    )
    resolved_relpaths_raw.update(alias_relpaths)
    resolved_relpaths = {
        str(key): Path(value) for key, value in resolved_relpaths_raw.items()
    }

    resolved_overrides: dict[str, Any] = dict(dict_file_overrides or {})
    resolved_overrides.update(alias_overrides)
    resolved_overrides = {key: value for key, value in resolved_overrides.items() if value}

    run_script_path = Path(run_script_relpath)
    normalized_pre_solve = tuple(pre_solve_commands or ())
    if _apply_case_mutation is None:
        from openfoam_driver.core.compatibility import legacy_generic_case_mutation

        _apply_case_mutation = legacy_generic_case_mutation

    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=OUTPUT_DIR_NAME,
    )

    normalized_cases = _normalize_case_specs(
        cases=cases,
        dict_file_overrides=resolved_overrides,
        dimension=dimension,
        parallel=parallel,
        touch_case_foam=touch_case_foam,
        openfoam_bashrc=openfoam_bashrc,
        solver_command=solver_command,
        pre_solve_commands=normalized_pre_solve,
    )

    # A case counts as generic when its *primary* declared dictionary file --
    # the first entry of ``dict_file_relpaths`` -- is absent from the folder.
    # Core imposes no vocabulary on that mapping: whichever file a caller (or
    # the legacy default in ``core.compatibility``) declares first is the one
    # whose presence marks the folder as belonging to that solver. Declaring no
    # dictionary files at all leaves the folder generic.
    primary_relpaths = list(resolved_relpaths.values())[:1]
    generic_case = (
        solver_command is None
        and cases is None
        and not resolved_overrides
        and not collect_patterns
        and not any((case_root / relpath).exists() for relpath in primary_relpaths)
    )

    return TutorialSpec(
        name=case_dir_name,
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=lambda: list(normalized_cases),
        apply_case=partial(
            _apply_case,
            dict_file_relpaths=resolved_relpaths,
            mutation_callback=_apply_case_mutation,
        ),
        metadata={
            "notes": "Core generic case runner for arbitrary tutorial folders.",
            "workflow_dag": _workflow_dag_for(
                solver_command=solver_command,
                pre_solve_commands=normalized_pre_solve,
            ),
            "dict_file_relpaths": {
                key: str(value) for key, value in resolved_relpaths.items()
            },
            "run_script_relpath": str(run_script_path),
            "collect_patterns": list(collect_patterns),
            "case_count": len(normalized_cases),
            "has_default_dict_file_overrides": bool(resolved_overrides),
            "solver_command": list(solver_command) if not isinstance(solver_command, str) and solver_command is not None else solver_command,
            "pre_solve_commands": list(pre_solve_commands or ()),
            "generic_case": generic_case,
        },
    )


def make_generic_case_spec(**kwargs: Any) -> TutorialSpec:
    """Solver-neutral alias used by registry case-folder discovery."""

    def _no_solver_mutation(*_args, **_kwargs) -> None:
        return None

    kwargs.setdefault("_apply_case_mutation", _no_solver_mutation)
    return make_spec(**kwargs)
