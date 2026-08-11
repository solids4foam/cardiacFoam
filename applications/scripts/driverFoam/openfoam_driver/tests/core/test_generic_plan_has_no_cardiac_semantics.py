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
from openfoam_driver.strict_planning import strict_plan

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
    """A plain OpenFOAM case with an Allrun and a workflow contract -- no
    cardiac dictionaries anywhere."""
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
    (case / "workflow_contract.json").write_text(
        json.dumps({"steps": [{"id": "run", "command": "Allrun", "depends_on": []}]})
    )
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
