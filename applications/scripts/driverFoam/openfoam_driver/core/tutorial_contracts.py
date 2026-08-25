#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     tutorial_contracts
#
# Description
#     Defines parameters and bounds for tutorial executions.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""CLI-only tutorial descriptors.

This module stays for CLI consumers (``introspection.py``, ``listVerifiers``,
``listIonicModels``). New automation should prefer Run documents; see
``schemas/run-document.json`` and the ``RunDocument`` model in
``openfoam_driver.core.runtime.run_model``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from .runtime.models import TutorialSpec

if TYPE_CHECKING:
    from .plugin_interface import DriverContext


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


def _find_tutorials_root(case_root: Path) -> Path:
    for candidate in (case_root.parent, *case_root.parents):
        if (candidate / "regressionTests").exists():
            return candidate
    return case_root.parent


def describe_tutorial_contract(
    spec: TutorialSpec,
    *,
    resolution: str,
    driver_context: "DriverContext",
) -> dict[str, Any]:
    case_root = spec.case_root
    tutorials_root = _find_tutorials_root(case_root)
    regression_root = tutorials_root / "regressionTests" / spec.name

    block_mesh_variants = _glob_relpaths(case_root / "system", "blockMeshDict*")
    reference_cases = []
    if regression_root.exists():
        reference_cases.append(
            str(regression_root.relative_to(tutorials_root))
        )

    # Split on the profile's own ``role``, not on a path prefix: the prefix
    # would make core re-derive plugin semantics from a string, and would
    # misfile a plugin-owned dictionary that happens to live under system/
    # (or a required initial-condition file that does not).
    required_rules = driver_context.capabilities.case_files.required_rules()
    core_required_files = tuple(
        rule.path for rule in required_rules
        if not rule.role.startswith("openfoam.")
    )
    solver_required_files = tuple(
        rule.path for rule in required_rules
        if rule.role.startswith("openfoam.")
    )
    conditional_files = tuple(
        driver_context.capabilities.case_files.conditional_files()
    )

    return {
        "name": spec.name,
        "resolution": resolution,
        "case_root": str(case_root),
        "setup_root": str(spec.setup_root),
        "output_dir": str(spec.output_dir),
        "core_required_files": _existing_relpaths(case_root, core_required_files),
        "solver_required_files": _existing_relpaths(case_root, solver_required_files),
        "conditional_files": _existing_relpaths(case_root, conditional_files),
        "mesh_files": block_mesh_variants,
        "constant_files": _glob_relpaths(case_root / "constant", "*"),
        "system_files": _glob_relpaths(case_root / "system", "*"),
        "reference_cases": reference_cases,
        "postprocess_modules": _glob_relpaths(case_root, "post_processing*.py"),
        "case_parameters": _case_parameter_contract(spec),
        "metadata": spec.metadata,
    }
