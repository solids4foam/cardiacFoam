from __future__ import annotations

import inspect
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

from .core.runtime.models import CaseConfig, TutorialSpec
from .core.runtime.registry import (
    list_available_tutorials,
    list_case_directories,
    list_tutorials,
    resolve_tutorial,
)
from .dict_entries import ELECTRO_PROPERTY_ENTRY_GROUPS, PHYSICS_PROPERTY_ENTRIES
from .gui_schema import describe_gui_schema
from .ionic_model_catalog import (
    ACTIVE_TENSION_MODEL_CATALOG,
    IONIC_MODEL_CATALOG,
    SOLVER_COMPATIBILITY_RULES,
)
from .launch import describe_launch_matrix
from .tutorial_contracts import describe_tutorial_contract

COMMON_OVERRIDE_KEYS = (
    "case_dir_name",
    "setup_dir_name",
    "output_dir_name",
    "run_script_relpath",
    "electro_property_overrides",
    "physics_property_overrides",
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
    if isinstance(value, set):
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


def _dict_entry_catalog() -> dict[str, Any]:
    return {
        "physicsProperties": [_serialize(asdict(entry)) for entry in PHYSICS_PROPERTY_ENTRIES],
        "electroProperties": {
            group_name: [_serialize(asdict(entry)) for entry in entries]
            for group_name, entries in ELECTRO_PROPERTY_ENTRY_GROUPS.items()
        },
    }


def _ionic_model_catalog() -> dict[str, Any]:
    return {
        "schema_version": "1.0",
        "ionic_models": {
            name: _serialize(asdict(entry))
            for name, entry in IONIC_MODEL_CATALOG.items()
        },
        "active_tension_models": {
            name: _serialize(asdict(entry))
            for name, entry in ACTIVE_TENSION_MODEL_CATALOG.items()
        },
        "solver_compatibility": [
            _serialize(rule) for rule in SOLVER_COMPATIBILITY_RULES
        ],
    }


def list_runs(runs_root: str | Path) -> list[dict[str, Any]]:
    import json
    root = Path(runs_root)
    manifests = []
    if not root.exists():
        return manifests

    for path in root.rglob("run_manifest.json"):
        try:
            payload = json.loads(path.read_text())
            manifests.append(payload)
        except Exception:
            pass

    manifests.sort(key=lambda m: m.get("started_at_utc") or "", reverse=True)
    return manifests


def describe_tutorial(
    tutorial: str,
    *,
    overrides: dict[str, Any] | None = None,
    config_path: str | Path | None = None,
    python_executable: str | None = None,
) -> dict[str, Any]:
    resolution = resolve_tutorial(tutorial, overrides=overrides)
    spec = resolution["factory"](**resolution["factory_overrides"])
    tutorials_root = Path(
        resolution["factory_overrides"].get("tutorials_root", spec.case_root.parent)
    )

    return {
        "requested_tutorial": tutorial,
        "resolution": resolution["resolution"],
        "resolved_name": resolution["resolved_name"],
        "registered_tutorials": list_tutorials(),
        "special_tutorial_aliases": list(SPECIAL_TUTORIAL_ALIASES),
        "available_tutorials": list_available_tutorials(tutorials_root),
        "case_directories": list_case_directories(tutorials_root),
        "common_override_keys": list(COMMON_OVERRIDE_KEYS),
        "make_spec": _describe_factory(resolution["factory"]),
        "factory_overrides": _serialize(resolution["factory_overrides"]),
        "spec": _describe_spec(spec),
        "tutorial_contract": _serialize(
            describe_tutorial_contract(
                spec,
                resolution=resolution["resolution"],
            )
        ),
        "dict_entries": _dict_entry_catalog(),
        "ionic_model_catalog": _ionic_model_catalog(),
        "gui_schema": describe_gui_schema(),
        "launch": describe_launch_matrix(
            tutorial,
            overrides=resolution["factory_overrides"],
            config_path=config_path,
            tutorials_root=tutorials_root,
            python_executable=python_executable,
        ),
    }
