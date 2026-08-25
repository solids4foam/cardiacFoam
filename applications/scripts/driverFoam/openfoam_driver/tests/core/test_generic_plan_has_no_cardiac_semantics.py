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
#     test_generic_plan_has_no_cardiac_semantics
#
# Description
#     Phase 1 exit gate: a plan produced under --plugin none must contain no
#     cardiacFoam command, field, required-file, utility, or override
#     semantics. Reading the code is not evidence -- this runs it.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import stat
from pathlib import Path

from openfoam_driver.core.plugin_interface import generic_openfoam_context
from openfoam_driver.core.introspection import describe_entry
from openfoam_driver.core.strict_planning import strict_plan

# Every token that would betray a cardiac assumption leaking into a plan
# produced for a non-cardiac solver.
_CARDIAC_TOKENS = (
    "cardiacFoam",
    "electroProperties",
    "physicsProperties",
    "ELECTRO_MODEL_COEFFS",
    "ionicModel",
    "myocardiumSolver",
    "activationTime",
    "phiE",
    "listCellModelsVariables",
    "bathBidomainInterfaceMetrics",
)


def _minimal_case(root: Path) -> Path:
    """A plain OpenFOAM case with an Allrun and no cardiac dictionaries."""
    case = root / "case"
    (case / "system").mkdir(parents=True)
    (case / "constant").mkdir(parents=True)
    (case / "system" / "controlDict").write_text(
        "FoamFile{version 2.0; format ascii; class dictionary; "
        "object controlDict;}\n"
        "application myGenericSolver;\nstartFrom startTime;\nstartTime 0;\n"
        "stopAt endTime;\nendTime 1;\ndeltaT 0.1;\nwriteControl timeStep;\n"
        "writeInterval 10;\n"
    )
    allrun = case / "Allrun"
    allrun.write_text("#!/bin/sh\necho generic-allrun-ran\n")
    allrun.chmod(allrun.stat().st_mode | stat.S_IEXEC)
    return case


def _generic_plan(tmp_path: Path) -> dict:
    case = _minimal_case(tmp_path)
    return strict_plan(
        str(case.relative_to(tmp_path)),
        overrides={"tutorials_root": str(tmp_path)},
        driver_context=generic_openfoam_context(),
    ).to_json()


def test_generic_plan_contains_no_cardiac_semantics(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("SKIP_ENV_DIAGNOSTICS", "1")
    monkeypatch.setenv("SKIP_MESH_DIAGNOSTICS", "1")
    blob = json.dumps(_generic_plan(tmp_path))
    leaked = [token for token in _CARDIAC_TOKENS if token in blob]
    assert leaked == [], f"cardiac semantics leaked into a generic plan: {leaked}"


def test_generic_plan_still_produces_a_usable_contract(tmp_path, monkeypatch) -> None:
    """Emptiness is not the goal -- the plan must still be runnable."""
    monkeypatch.setenv("SKIP_ENV_DIAGNOSTICS", "1")
    monkeypatch.setenv("SKIP_MESH_DIAGNOSTICS", "1")
    payload = _generic_plan(tmp_path)
    assert payload["workflow_dag"]["steps"], "generic plan must have runnable steps"
    assert payload["capability_manifest"]["allowed_commands"]["utilities"] == {}
    assert "cardiacFoam" not in payload["capability_manifest"]["allowed_commands"]["core"]


def test_generic_describe_override_surface_has_no_cardiac_semantics(
    tmp_path, monkeypatch
) -> None:
    """The spec's exit gate names "override semantics", but those live in the
    describe payload (``config_schema``, ``dict_entries``) -- ``strict_plan``
    does not emit them, so gating only on the plan left the one clause naming
    the override surface checked against a payload that cannot contain it."""
    monkeypatch.setenv("SKIP_ENV_DIAGNOSTICS", "1")
    monkeypatch.setenv("SKIP_MESH_DIAGNOSTICS", "1")
    case = _minimal_case(tmp_path)
    payload = describe_entry(
        str(case.relative_to(tmp_path)),
        overrides={"tutorials_root": str(tmp_path)},
        driver_context=generic_openfoam_context(),
    )
    override_surface = {
        "config_schema": payload["config_schema"],
        "dict_entries": payload["dict_entries"],
    }
    blob = json.dumps(override_surface)
    leaked = [token for token in _CARDIAC_TOKENS if token in blob]
    assert leaked == [], f"cardiac override semantics leaked: {leaked}"


def test_generic_spec_metadata_names_dict_files_generically(
    tmp_path, monkeypatch
) -> None:
    """Phase 2 (P2.6) closes the residual that Phase 1 only documented.

    Core's generic-case metadata used to carry two hard-coded cardiac field
    names (``electro_properties_relpath``/``physics_properties_relpath``) even
    under ``--plugin none``. It now carries a single generic
    ``dict_file_relpaths`` mapping whose *keys* are chosen by whoever declares
    the dictionaries; the cardiac-shaped default values still arrive via the
    named ``core.compatibility`` seam, so the paths are unchanged while the
    core-owned vocabulary is not."""
    monkeypatch.setenv("SKIP_ENV_DIAGNOSTICS", "1")
    monkeypatch.setenv("SKIP_MESH_DIAGNOSTICS", "1")
    case = _minimal_case(tmp_path)
    payload = describe_entry(
        str(case.relative_to(tmp_path)),
        overrides={"tutorials_root": str(tmp_path)},
        driver_context=generic_openfoam_context(),
    )
    metadata = payload["spec"]["metadata"]
    assert "electro_properties_relpath" not in metadata
    assert "physics_properties_relpath" not in metadata
    assert "has_default_electro_property_overrides" not in metadata
    assert "has_default_physics_property_overrides" not in metadata
    assert metadata["has_default_dict_file_overrides"] is False
    assert set(metadata["dict_file_relpaths"]) == {"electro", "physics"}
