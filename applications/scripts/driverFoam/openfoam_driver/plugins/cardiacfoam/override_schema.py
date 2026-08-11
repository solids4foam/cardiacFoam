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
#     override_schema
#
# Description
#     cardiacFoam's authored configuration vocabulary: the --config schema
#     prose and worked examples, and the electro/physics document shape of the
#     dictionary-entry catalog. Core assembles and serializes these; the
#     vocabulary itself is solver knowledge and lives here.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from typing import Any

# Spec keys that describe where a case lives rather than what it computes.
# Excluded from the worked example so the schema advertises physics knobs.
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


def dict_entry_catalog(catalog: Any, groups: dict[str, Any]) -> dict[str, Any]:
    """Return the cardiac document shape of the dictionary-entry catalog.

    ``physicsProperties`` is a flat sequence while ``electroProperties`` is
    grouped; that asymmetry mirrors the two OpenFOAM dictionaries this solver
    reads and is deliberately not a core convention. Values are returned
    unserialized -- core owns serialization.
    """
    return {
        "physicsProperties": list(catalog.entries_for("physicsProperties")),
        "electroProperties": {
            group_name: list(entries) for group_name, entries in groups.items()
        },
    }


def config_schema(
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
