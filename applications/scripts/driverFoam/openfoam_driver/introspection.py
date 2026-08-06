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
from typing import Any

def __get_capabilities():
    from openfoam_driver.core.plugin_interface import get_active_plugin
    return get_active_plugin().get_capabilities()


from .core.runtime.models import CaseConfig, TutorialSpec
from .core.runtime.registry import (
    list_entries,
    list_available_tutorials,
    list_case_directories,
    list_tutorials,
    resolve_entry,
)
from .capability_manifest import build_capability_manifest, resolve_case_models
from .core.runtime.execution_context import resolve_execution_context
from .dict_entries import get_electro_property_entry_groups, PHYSICS_PROPERTY_ENTRIES
from .strict_planning import _run_launch_description
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
            for group_name, entries in get_electro_property_entry_groups().items()
        },
    }


def _ionic_model_catalog() -> dict[str, Any]:
    return {
        "schema_version": "1.0",
        "ionic_models": {
            name: _serialize(asdict(entry))
            for name, entry in __get_capabilities().get("ionic_models", {}).items()
        },
        "solver_compatibility": [
            _serialize(rule) for rule in __get_capabilities().get("solver_compatibility_rules", [])
        ],
    }


def _active_tension_catalog() -> dict[str, Any]:
    return {
        "schema_version": "1.0",
        "active_tension_models": {
            name: _serialize(asdict(entry))
            for name, entry in __get_capabilities().get("active_tension_models", {}).items()
        },
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
            "or wrapped in an entry-named key (multi-entry files where different "
            "sections apply to different entries)."
        ),
        "top_level_shapes": {
            "flat": {
                "description": (
                    "A single JSON object whose keys are spec parameters and/or "
                    "common override keys. Applies to the entry named on --entry."
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
                    "A JSON object keyed by entry name. Each value is a flat "
                    "section. Use this when one file covers multiple entries. "
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
                    "Parameters accepted by make_spec() for this entry. "
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
    """Static schema description for run_manifest.json.

    This describes the legacy sim/post/all CLI's own manifest -- see
    payload["strict_launch"] for the canonical way to actually execute an
    entry today (run --strict), whose own state lives in
    output_dir/workflow_state.json instead, written by every workflow step
    as it runs."""
    return {
        "description": (
            "run_manifest.json is the run-state source of truth. "
            "It is written to output_dir/run_manifest.json and updated after every "
            "case completes. Poll this file to track run progress."
        ),
        "schema_version": "2.3",
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
            "schema_version": "string — manifest format version (currently '2.3'; v2.x is additive-only)",
            "run_id": "string — unique ID for this run (timestamp + random suffix)",
            "requested_action": "string — 'sim', 'post', or 'all'",
            "entry": "string — selected entry name",
            "entry_kind": "string | null — entry classification such as registered_tutorial or workflow_case",
            "entry_path": "string | null — relative path of the resolved entry under tutorials/",
            "source_type": "string | null — spec_factory, workflow_contract, workflow_reference_case, filesystem_case, or generic_alias",
            "workflow_family": "string | null — workflow family name when the entry belongs to one",
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
            "artifacts_manifest_path": "string | null — path to artifacts_manifest.json (predicted DataArtifacts for the current case state; v2.2+)",
            "artifacts_realized_path": "string | null — path to artifacts_realized.json (v1.1: cases[] array, one entry per sweep case; predicted-vs-actual reconciliation; written only at terminal status on non-dry runs; v2.3+)",
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
        family = families.setdefault(
            family_name,
            {
                "workflow_family": family_name,
                "template_entry": None,
                "reference_cases": [],
                "workflow_templates": [],
            },
        )

        entry_kind = str(entry["entry_kind"])
        if entry_kind == "workflow_template":
            authoring_contract = (
                Path(tutorials_root) / str(entry["entry_path"]) / "workflow_contract.json"
            )
            payload = None
            if authoring_contract.exists():
                import json

                payload = json.loads(authoring_contract.read_text())
            family["template_entry"] = {
                "entry_name": entry["entry_name"],
                "entry_path": entry["entry_path"],
                "entry_kind": entry_kind,
                "is_runnable": entry["is_runnable"],
            }
            family["workflow_templates"] = list((payload or {}).get("workflow_templates", []))
        elif entry_kind == "workflow_case":
            family["reference_cases"].append(
                {
                    "entry_name": entry["entry_name"],
                    "entry_path": entry["entry_path"],
                    "entry_kind": entry_kind,
                    "is_runnable": entry["is_runnable"],
                }
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
) -> dict[str, Any]:
    resolution = resolve_entry(entry, entry_kind=entry_kind, overrides=overrides)
    spec = resolution["factory"](**resolution["factory_overrides"])
    tutorials_root = Path(
        resolution["factory_overrides"].get("tutorials_root", spec.case_root.parent)
    )
    entry_catalog = list_entries(tutorials_root)
    workflow_catalog = _workflow_catalog(tutorials_root, entry_catalog)

    make_spec_info = _describe_factory(resolution["factory"])
    _solver, _ionic, _active_tension = resolve_case_models(spec.case_root)
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
        "active_tension_catalog": _active_tension_catalog(),
        "strict_launch": _run_launch_description(
            resolution["resolved_name"],
            resolve_execution_context(spec),
            entry_kind=resolution["entry_kind"],
            config_path=config_path,
        ),
        "config_schema": _describe_config_schema(
            resolution["resolved_name"],
            make_spec_info,
        ),
        "manifest_schema": _manifest_schema(),
        "capability_manifest": build_capability_manifest(
            resolved_solver=_solver,
            resolved_ionic_model=_ionic,
            resolved_active_tension=_active_tension,
        ),
    }


def describe_tutorial(
    tutorial: str,
    *,
    overrides: dict[str, Any] | None = None,
    config_path: str | Path | None = None,
) -> dict[str, Any]:
    return describe_entry(
        tutorial,
        overrides=overrides,
        config_path=config_path,
    )
