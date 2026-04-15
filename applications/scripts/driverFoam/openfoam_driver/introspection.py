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


_INFRASTRUCTURE_SPEC_KEYS = frozenset(
    {
        "tutorials_root",
        "case_dir_name",
        "setup_dir_name",
        "output_dir_name",
        "run_script_relpath",
        "postprocess_strict_artifacts",
    }
)


def _describe_config_schema(
    tutorial_name: str,
    make_spec_info: dict[str, Any],
) -> dict[str, Any]:
    """Build the config_schema payload for a tutorial.

    Returns a machine-readable description of the --config JSON format accepted
    by the driver, including a worked example specific to this tutorial.
    """
    # Collect spec-level parameters (exclude infrastructure keys)
    spec_params = {
        k: v
        for k, v in make_spec_info.get("parameters", {}).items()
        if k not in _INFRASTRUCTURE_SPEC_KEYS
    }

    # Build a minimal worked example
    example_section: dict[str, Any] = {}
    # Pick the first non-required spec param that has a readable default for demo
    for param_name, param_info in spec_params.items():
        if "default" in param_info and param_info["default"] is not None:
            example_section[param_name] = param_info["default"]
            break
    # Always show an electro_property_overrides example with real driver_path keys
    example_section["electro_property_overrides"] = {
        "$ELECTRO_MODEL_COEFFS.initialODEStep": "1e-5",
        "$ELECTRO_MODEL_COEFFS.maxSteps": "1000",
    }

    return {
        "description": (
            "Describes the --config JSON file format accepted by openfoam_driver. "
            "The config is a JSON object. It may be flat (applies to one tutorial) "
            "or wrapped in a tutorial-named key (multi-tutorial files where different "
            "sections apply to different tutorials)."
        ),
        "top_level_shapes": {
            "flat": {
                "description": (
                    "A single JSON object whose keys are spec parameters and/or "
                    "common override keys. Applies to the tutorial named on --tutorial."
                ),
                "example_snippet": {
                    "ionic_models": ["TNNP"],
                    "electro_property_overrides": {
                        "$ELECTRO_MODEL_COEFFS.chi": "140000",
                    },
                },
            },
            "wrapped": {
                "description": (
                    "A JSON object keyed by tutorial name. Each value is a flat "
                    "section. Use this when one file covers multiple tutorials. "
                    "Key matching is case-insensitive."
                ),
                "example_snippet": {
                    tutorial_name: {
                        "ionic_models": ["TNNP"],
                        "electro_property_overrides": {
                            "$ELECTRO_MODEL_COEFFS.chi": "140000",
                        },
                    }
                },
            },
        },
        "section_fields": {
            "spec_parameters": {
                "description": (
                    "Parameters accepted by make_spec() for this tutorial. "
                    "These are the high-level knobs (e.g. ionic_models, n_beats, "
                    "solvers). Place them at the top level of the config section."
                ),
                "available_keys": list(spec_params.keys()),
            },
            "electro_property_overrides": {
                "description": (
                    "Overrides for entries in constant/electroProperties. "
                    "Keys are driver_path strings from dict_entries.electroProperties. "
                    "The $ELECTRO_MODEL_COEFFS token is resolved automatically to "
                    "the actual solver coeffs dict (e.g. monodomainSolverCoeffs)."
                ),
                "shorthand_format": {
                    "description": (
                        "Recommended. A mapping from driver_path to value string. "
                        "Use exactly the driver_path values listed in dict_entries."
                    ),
                    "example": {
                        "$ELECTRO_MODEL_COEFFS.chi": "140000",
                        "$ELECTRO_MODEL_COEFFS.cm": "0.01",
                        "$ELECTRO_MODEL_COEFFS.initialODEStep": "1e-5",
                        "$ELECTRO_MODEL_COEFFS.ionicModel": "TNNP",
                    },
                },
                "explicit_format": {
                    "description": (
                        "A list of {key, scope, value} objects. Use when you need "
                        "to address a specific sub-dictionary by name without the "
                        "$ELECTRO_MODEL_COEFFS shorthand."
                    ),
                    "example": [
                        {
                            "key": "chi",
                            "scope": ["monodomainSolverCoeffs"],
                            "value": "140000",
                        },
                        {
                            "key": "ionicModel",
                            "scope": ["monodomainSolverCoeffs"],
                            "value": "TNNP",
                        },
                    ],
                },
                "note": (
                    "A single config section may supply either format but not both "
                    "simultaneously. The shorthand form is sufficient for all "
                    "driver_path entries listed in dict_entries."
                ),
            },
            "physics_property_overrides": {
                "description": (
                    "Overrides for constant/physicsProperties. Same format as "
                    "electro_property_overrides (shorthand mapping or explicit list). "
                    "No $ELECTRO_MODEL_COEFFS token — keys are bare property names."
                ),
                "example": {"type": "electroMechanicalModel"},
            },
            "case_dir_name": {
                "type": "string",
                "description": "Override the tutorial case directory name.",
            },
            "setup_dir_name": {
                "type": "string",
                "description": (
                    "Override the setup directory. Defaults to "
                    "<case_dir_name>_setup if omitted."
                ),
            },
            "output_dir_name": {
                "type": "string",
                "description": "Override where run outputs and the manifest are written.",
            },
            "run_script_relpath": {
                "type": "string",
                "description": "Relative path to a different run script within the case.",
            },
        },
        "worked_example": {
            "description": f"Minimal working config for the '{tutorial_name}' tutorial.",
            "json": {tutorial_name: example_section},
        },
    }


def _manifest_schema() -> dict[str, Any]:
    """Static schema description for run_manifest.json."""
    return {
        "description": (
            "run_manifest.json is the run-state source of truth. "
            "It is written to output_dir/run_manifest.json and updated after every "
            "case completes. Poll this file to track run progress."
        ),
        "schema_version": "2.0",
        "file_location": "output_dir/run_manifest.json  (see launch.<action>.manifest_path)",
        "companion_file": (
            "output_dir/action_events.jsonl — append-only JSONL log with one "
            "event per line: sim_started, case_started, case_finished, "
            "sim_finished, postprocess_started, postprocess_finished, all_started, all_finished."
        ),
        "polling_guidance": (
            "Poll every 15-30 seconds. Stop when status is one of the terminal states. "
            "Reading the file is safe at any time — it is written atomically."
        ),
        "top_level_fields": {
            "schema_version": "string — manifest format version (currently '2.0')",
            "run_id": "string — unique ID for this run (timestamp + random suffix)",
            "requested_action": "string — 'sim', 'post', or 'all'",
            "tutorial": "string — tutorial name",
            "case_root": "string — absolute path to the case directory",
            "setup_root": "string — absolute path to the setup directory",
            "output_dir": "string — absolute path to the output directory",
            "dry_run": "boolean",
            "continue_on_error": "boolean",
            "status": "string — see status_values below",
            "postprocess_status": "string — see postprocess_status_values below",
            "current_case_id": "string | null — case currently executing, null between cases",
            "started_at_utc": "string | null — ISO 8601 UTC timestamp",
            "updated_at_utc": "string — ISO 8601 UTC timestamp of last write",
            "finished_at_utc": "string | null — ISO 8601 UTC timestamp, null until terminal",
            "total_cases": "integer",
            "planned_cases": "integer — cases with status 'planned' (dry_run only)",
            "completed_cases": "integer — cases with status 'ok'",
            "failed_cases": "integer — cases with status 'failed'",
            "error": "string | null — top-level error message if run failed early",
            "plots_manifest_path": "string | null — path to plots.json if postprocess produced plots",
            "human_report_path": "string — path to run_report.md",
            "results": "array of CaseResult objects — see case_result_fields",
        },
        "status_values": {
            "running": "Simulation is in progress.",
            "completed": "All cases finished successfully. Terminal.",
            "completed_with_failures": "All cases ran; at least one failed. Terminal.",
            "failed": "A case failed and continue_on_error=false. Terminal.",
            "postprocessing": "Simulations done; postprocessing is now running.",
            "postprocess_failed": "Postprocessing raised an exception. Terminal.",
            "planned": "Dry-run completed — no actual simulation was run. Terminal.",
        },
        "postprocess_status_values": {
            "not_started": "Postprocessing has not begun.",
            "running": "Postprocessing is executing.",
            "completed": "Postprocessing finished successfully.",
            "failed": "Postprocessing raised an exception.",
            "skipped": "Dry-run mode; postprocessing was not attempted.",
        },
        "terminal_states": [
            "completed",
            "completed_with_failures",
            "failed",
            "postprocess_failed",
            "planned",
        ],
        "case_result_fields": {
            "case_id": "string — unique case identifier",
            "status": "'ok' | 'failed' | 'planned'",
            "duration_s": "float — wall-clock seconds",
            "params": "object — parameter values for this case",
            "error": "string | null — exception message if status='failed'",
            "index": "integer — 1-based position in the case list",
            "total_cases": "integer",
            "started_at_utc": "string | null — ISO 8601",
            "finished_at_utc": "string | null — ISO 8601",
        },
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

    make_spec_info = _describe_factory(resolution["factory"])
    return {
        "requested_tutorial": tutorial,
        "resolution": resolution["resolution"],
        "resolved_name": resolution["resolved_name"],
        "registered_tutorials": list_tutorials(),
        "special_tutorial_aliases": list(SPECIAL_TUTORIAL_ALIASES),
        "available_tutorials": list_available_tutorials(tutorials_root),
        "case_directories": list_case_directories(tutorials_root),
        "common_override_keys": list(COMMON_OVERRIDE_KEYS),
        "make_spec": make_spec_info,
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
        "config_schema": _describe_config_schema(
            resolution["resolved_name"],
            make_spec_info,
        ),
        "manifest_schema": _manifest_schema(),
    }
