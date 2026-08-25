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
#     introspection
#
# Description
#     Provides reflection tools to query runtime configurations.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import inspect
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .plugin_interface import DriverContext


from .runtime.models import CaseConfig, TutorialSpec
from .runtime.registry import (
    list_entries,
    list_available_tutorials,
    list_case_directories,
    list_tutorials,
    resolve_entry,
)
from .runtime.execution_context import resolve_execution_context
from openfoam_driver.core.strict_planning import _run_launch_description
from openfoam_driver.core.tutorial_contracts import describe_tutorial_contract

COMMON_OVERRIDE_KEYS = (
    "case_dir_name",
    "setup_dir_name",
    "output_dir_name",
    "run_script_relpath",
    "dict_file_relpaths",
    "dict_file_overrides",
    "postprocess_strict_artifacts",
)

SPECIAL_TUTORIAL_ALIASES = ("genericCase", "randomCase")


def _serialize(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return _serialize(asdict(value))
    if isinstance(value, dict):
        return {str(key): _serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize(item) for item in value]
    if isinstance(value, (set, frozenset)):
        # frozenset is NOT a subclass of set; without it, DictEntry.phases fell
        # through to repr() and shipped "frozenset({'physics'})" as JSON.
        return sorted(_serialize(item) for item in value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def _annotation_to_string(annotation: Any) -> str | None:
    if annotation is inspect.Signature.empty:
        return None
    if isinstance(annotation, str):
        return annotation
    return repr(annotation).replace("typing.", "")


def _describe_parameter(parameter: inspect.Parameter) -> dict[str, Any]:
    payload = {
        "kind": parameter.kind.name.lower(),
        "required": parameter.default is inspect.Signature.empty,
    }
    annotation = _annotation_to_string(parameter.annotation)
    if annotation is not None:
        payload["annotation"] = annotation
    if parameter.default is not inspect.Signature.empty:
        payload["default"] = _serialize(parameter.default)
    return payload


def _describe_factory(factory: object) -> dict[str, Any]:
    if not callable(factory):
        raise TypeError(f"Factory is not callable: {factory!r}")
    signature = inspect.signature(factory)
    return {
        "callable": f"{factory.__module__}.{factory.__name__}",
        "parameters": {
            name: _describe_parameter(parameter)
            for name, parameter in signature.parameters.items()
        },
    }


def _describe_cases(cases: list[CaseConfig]) -> dict[str, Any]:
    return {
        "count": len(cases),
        "items": [
            {
                "case_id": case.case_id,
                "params": _serialize(case.params),
            }
            for case in cases
        ],
    }


def _describe_spec(spec: TutorialSpec) -> dict[str, Any]:
    cases = spec.build_cases()
    return {
        "name": spec.name,
        "case_root": str(spec.case_root),
        "setup_root": str(spec.setup_root),
        "output_dir": str(spec.output_dir),
        "metadata": _serialize(spec.metadata),
        "cases": _describe_cases(cases),
    }


def _dict_entry_catalog(driver_context: "DriverContext | None" = None) -> dict[str, Any]:
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    # The document names and their shape are plugin vocabulary; core only
    # serializes whatever structure the plugin declares.
    return _serialize(
        driver_context.capabilities.override_schema.dict_entry_catalog()
    )


def _plugin_catalogs(driver_context: "DriverContext | None" = None) -> dict[str, Any]:
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    # The catalog names and their contents are plugin vocabulary (e.g. the
    # cardiac plugin's ionic_model_catalog/active_tension_catalog); core only
    # namespaces the whole mapping under this key and serializes it.
    return _serialize(
        dict(driver_context.capabilities.named_catalogs.catalogs())
    )


def _describe_config_schema(
    tutorial_name: str,
    make_spec_info: dict[str, Any],
    driver_context: "DriverContext | None" = None,
) -> dict[str, Any]:
    """Return the plugin-authored config_schema payload for a tutorial.

    The vocabulary (override tokens, examples, document names) is solver
    knowledge and lives in the active plugin; core only routes the request.
    """
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    return driver_context.capabilities.override_schema.config_schema(
        tutorial_name, make_spec_info,
    )


def _run_state_schema() -> dict[str, Any]:
    """Static schema description for workflow_state.json."""
    return {
        "description": (
            "workflow_state.json is the run-state source of truth. It is "
            "written to output_dir/workflow_state.json by the strict workflow "
            "orchestrator and updated after every workflow step. Poll this "
            "file to track run progress."
        ),
        "schema_version": "3.0",
        "file_location": (
            "output_dir/workflow_state.json  (see launch.<action>.workflowStatePath)"
        ),
        "polling_guidance": (
            "Poll every 15-30 seconds. Read status and current_step_id; each "
            "step carries its own status, attempts, exit code, and log path. "
            "Stop when status is a terminal state."
        ),
    }


def _workflow_catalog(
    tutorials_root: Path,
    entry_catalog: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    families: dict[str, dict[str, Any]] = {}
    for entry in entry_catalog:
        workflow_family = entry.get("workflow_family")
        if not workflow_family:
            continue
        family_name = str(workflow_family)
        families.setdefault(
            family_name,
            {
                "workflow_family": family_name,
                "template_entry": None,
                "reference_cases": [],
                "workflow_templates": [],
            },
        )

    return sorted(families.values(), key=lambda item: item["workflow_family"].casefold())


def _matching_workflow(
    workflow_catalog: list[dict[str, Any]],
    workflow_family: str | None,
) -> dict[str, Any] | None:
    if not workflow_family:
        return None
    for family in workflow_catalog:
        if family["workflow_family"] == workflow_family:
            return family
    return None


def describe_entry(
    entry: str,
    *,
    entry_kind: str | None = None,
    overrides: dict[str, Any] | None = None,
    config_path: str | Path | None = None,
    driver_context: "DriverContext | None" = None,
) -> dict[str, Any]:
    from .compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    resolution = resolve_entry(
        entry,
        entry_kind=entry_kind,
        overrides=overrides,
        driver_context=driver_context,
    )
    spec = resolution["factory"](**resolution["factory_overrides"])
    tutorials_root = Path(
        resolution["factory_overrides"].get("tutorials_root", spec.case_root.parent)
    )
    entry_catalog = list_entries(tutorials_root, driver_context=driver_context)
    workflow_catalog = _workflow_catalog(tutorials_root, entry_catalog)

    make_spec_info = _describe_factory(resolution["factory"])
    return {
        "requested_entry": entry,
        "resolution": resolution["resolution"],
        "resolved_name": resolution["resolved_name"],
        "entry": {
            "entry_name": resolution["entry_name"],
            "entry_kind": resolution["entry_kind"],
            "entry_path": resolution["entry_path"],
            "is_runnable": resolution["is_runnable"],
            "source_type": resolution["source_type"],
            "workflow_family": resolution["workflow_family"],
        },
        "entry_kind": resolution["entry_kind"],
        "entry_catalog": _serialize(entry_catalog),
        "workflow": _serialize(
            _matching_workflow(workflow_catalog, resolution["workflow_family"])
        ),
        "workflow_catalog": _serialize(workflow_catalog),
        "is_runnable": resolution["is_runnable"],
        "registered_tutorials": list_tutorials(driver_context),
        "special_tutorial_aliases": list(SPECIAL_TUTORIAL_ALIASES),
        "available_tutorials": list_available_tutorials(
            tutorials_root, driver_context=driver_context,
        ),
        "case_directories": list_case_directories(
            tutorials_root, driver_context=driver_context,
        ),
        "common_override_keys": list(COMMON_OVERRIDE_KEYS),
        "make_spec": make_spec_info,
        "factory_overrides": _serialize(resolution["factory_overrides"]),
        "spec": _describe_spec(spec),
        "tutorial_contract": _serialize(
            describe_tutorial_contract(
                spec,
                resolution=resolution["resolution"],
                driver_context=driver_context,
            )
        ),
        "dict_entries": _dict_entry_catalog(driver_context),
        "plugin_catalogs": _plugin_catalogs(driver_context),
        "strict_launch": _run_launch_description(
            resolution["resolved_name"],
            resolve_execution_context(spec),
            entry_kind=resolution["entry_kind"],
            config_path=config_path,
        ),
        "config_schema": _describe_config_schema(
            resolution["resolved_name"],
            make_spec_info,
            driver_context,
        ),
        "run_state_schema": _run_state_schema(),
        "capability_manifest": _serialize({
            **dict(driver_context.capabilities.manifest.manifest()),
            "plugin_identity": driver_context.identity.to_json(),
        }),
    }


def describe_tutorial(
    tutorial: str,
    *,
    overrides: dict[str, Any] | None = None,
    config_path: str | Path | None = None,
    driver_context: "DriverContext | None" = None,
) -> dict[str, Any]:
    return describe_entry(
        tutorial,
        overrides=overrides,
        config_path=config_path,
        driver_context=driver_context,
    )


def describe_launch_matrix(
    driver_context: "DriverContext | None" = None,
) -> list[dict[str, Any]]:
    """Return every registered entry. Alias kept for AGENT_GUIDE.md compatibility.

    Equivalent to list_entries() from openfoam_driver.core.runtime.registry.
    """
    return list_entries(driver_context=driver_context)
