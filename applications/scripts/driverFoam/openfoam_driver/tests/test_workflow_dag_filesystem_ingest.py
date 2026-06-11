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
#     test_workflow_dag_filesystem_ingest
#
# Description
#     Tests workflow dag filesystem ingest logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for workflow_dag ingest from on-disk workflow_contract.json (plan §6.2).

Filesystem cases (loaded via registry.py) must have their workflow_dag populated
from the on-disk workflow_contract.json when that file contains a 'steps' array.
"""
from __future__ import annotations

import json
import unittest
from pathlib import Path
import tempfile

from openfoam_driver.core.runtime.registry import load_tutorial_spec


class TestWorkflowDagFilesystemIngest(unittest.TestCase):
    """workflow_dag is populated from workflow_contract.json for filesystem cases."""

    def _write_case_files(self, case_root: Path) -> None:
        (case_root / "constant").mkdir(parents=True, exist_ok=True)
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver singleCellSolver;\n"
        )
        (case_root / "constant" / "physicsProperties").write_text(
            "type electroModel;\n"
        )

    def test_filesystem_case_workflow_dag_matches_workflow_contract(self) -> None:
        """spec.metadata['workflow_dag'] round-trips the steps from workflow_contract.json."""
        steps = [
            {"id": "mesh", "command": "blockMesh", "depends_on": []},
            {"id": "solve", "command": "cardiacFoam", "depends_on": ["mesh"]},
        ]
        contract = {
            "tutorial_family": "custom",
            "status": {"runnable_without_substitution": True},
            "steps": steps,
        }

        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "myCase"
            self._write_case_files(case_root)
            (case_root / "workflow_contract.json").write_text(json.dumps(contract))

            spec = load_tutorial_spec(
                "myCase",
                overrides={"tutorials_root": tutorials_root},
            )

            dag = spec.metadata.get("workflow_dag")
            self.assertIsNotNone(dag, "workflow_dag must be populated from workflow_contract.json")
            self.assertIn("steps", dag)
            self.assertEqual(dag["steps"], steps)

    def test_filesystem_case_without_contract_has_no_workflow_dag(self) -> None:
        """Filesystem case with no workflow_contract.json → workflow_dag absent or None."""
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "bareCase"
            self._write_case_files(case_root)
            # No workflow_contract.json written

            spec = load_tutorial_spec(
                "bareCase",
                overrides={"tutorials_root": tutorials_root},
            )

            # Either missing or None — both are acceptable for the None path
            dag = spec.metadata.get("workflow_dag")
            self.assertIsNone(dag, "workflow_dag must be None when workflow_contract.json is absent")

    def test_filesystem_case_contract_without_steps_has_no_workflow_dag(self) -> None:
        """workflow_contract.json without a 'steps' key → workflow_dag is None."""
        contract = {
            "tutorial_family": "custom",
            "status": {"runnable_without_substitution": True},
        }

        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "noStepsCase"
            self._write_case_files(case_root)
            (case_root / "workflow_contract.json").write_text(json.dumps(contract))

            spec = load_tutorial_spec(
                "noStepsCase",
                overrides={"tutorials_root": tutorials_root},
            )

            dag = spec.metadata.get("workflow_dag")
            self.assertIsNone(dag, "workflow_dag must be None when contract has no 'steps' key")


if __name__ == "__main__":
    unittest.main()
