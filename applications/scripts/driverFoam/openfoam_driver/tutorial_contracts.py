from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .core.runtime.models import TutorialSpec


CORE_REQUIRED_FILES = (
    "constant/electroProperties",
    "constant/physicsProperties",
)

SOLVER_REQUIRED_FILES = (
    "system/controlDict",
    "system/fvSchemes",
    "system/fvSolution",
)

CONDITIONAL_FILES = (
    "system/decomposeParDict",
    "system/blockMeshDict",
    "Allrun",
    "Allclean",
    "README.md",
    "runRegressionTest.sh",
)


def _existing_relpaths(case_root: Path, candidates: tuple[str, ...]) -> list[str]:
    existing: list[str] = []
    for relpath in candidates:
        if (case_root / relpath).exists():
            existing.append(relpath)
    return existing


def _glob_relpaths(root: Path, pattern: str) -> list[str]:
    if not root.exists():
        return []
    return sorted(
        str(path.relative_to(root))
        for path in root.glob(pattern)
        if not path.name.startswith(".")
    )


def _unique_case_param_values(spec: TutorialSpec, key: str) -> list[Any]:
    values = []
    seen: set[str] = set()
    for case in spec.build_cases():
        if key not in case.params:
            continue
        value = case.params[key]
        marker = repr(value)
        if marker in seen:
            continue
        seen.add(marker)
        values.append(value)
    return values


def _case_parameter_contract(spec: TutorialSpec) -> dict[str, list[Any]]:
    cases = spec.build_cases()
    if not cases:
        return {}

    keys = sorted({key for case in cases for key in case.params})
    return {
        key: _unique_case_param_values(spec, key)
        for key in keys
    }


def _read_json_if_exists(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _find_tutorials_root(case_root: Path) -> Path:
    for candidate in (case_root.parent, *case_root.parents):
        if (candidate / "regressionTests").exists():
            return candidate
    return case_root.parent


def describe_tutorial_contract(
    spec: TutorialSpec,
    *,
    resolution: str,
) -> dict[str, Any]:
    case_root = spec.case_root
    tutorials_root = _find_tutorials_root(case_root)
    regression_root = tutorials_root / "regressionTests" / spec.name

    block_mesh_variants = _glob_relpaths(case_root / "system", "blockMeshDict*")
    authoring_contract_path = case_root / "workflow_contract.json"
    authoring_contract = _read_json_if_exists(authoring_contract_path)
    reference_cases = []
    if regression_root.exists():
        reference_cases.append(
            str(regression_root.relative_to(tutorials_root))
        )

    return {
        "name": spec.name,
        "resolution": resolution,
        "case_root": str(case_root),
        "setup_root": str(spec.setup_root),
        "output_dir": str(spec.output_dir),
        "core_required_files": _existing_relpaths(case_root, CORE_REQUIRED_FILES),
        "solver_required_files": _existing_relpaths(case_root, SOLVER_REQUIRED_FILES),
        "conditional_files": _existing_relpaths(case_root, CONDITIONAL_FILES),
        "mesh_files": block_mesh_variants,
        "constant_files": _glob_relpaths(case_root / "constant", "*"),
        "system_files": _glob_relpaths(case_root / "system", "*"),
        "reference_cases": reference_cases,
        "postprocess_modules": _glob_relpaths(case_root, "post_processing*.py"),
        "case_parameters": _case_parameter_contract(spec),
        "authoring_contract_path": (
            str(authoring_contract_path.relative_to(case_root))
            if authoring_contract_path.exists()
            else None
        ),
        "authoring_contract": authoring_contract,
        "metadata": spec.metadata,
    }
