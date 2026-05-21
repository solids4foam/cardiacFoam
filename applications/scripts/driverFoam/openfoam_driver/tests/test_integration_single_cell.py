"""End-to-end integration test against a real cardiacFoam binary.

This test is marked ``@pytest.mark.integration`` — it is skipped unless
the cardiacFoam binary is available on PATH. Run explicitly with::

    pytest -m integration openfoam_driver/tests/test_integration_single_cell.py -v

The test exercises the full agent pipeline:

  build_and_launch(...)
    -> build_electro_properties + build_physics_properties
    -> write to disk
    -> DriverEngine.run_simulations()
    -> cardiacFoam binary executes
    -> engine writes run_manifest, artifacts_manifest, artifacts_realized

It is the proof-of-life test for the autonomous-agent contract.
"""
from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pytest


_CARDIACFOAM_AVAILABLE = shutil.which("cardiacFoam") is not None


@pytest.mark.integration
@unittest.skipUnless(_CARDIACFOAM_AVAILABLE, "cardiacFoam binary not on PATH")
class TestSingleCellEndToEnd(unittest.TestCase):
    """Build a singleCell AlievPanfilov case from intent, run it, and
    verify the realized manifest reports at least one matched artifact."""

    def test_build_and_launch_produces_matched_artifacts(self) -> None:
        from openfoam_driver.specs.dict_builder import build_and_launch

        with tempfile.TemporaryDirectory() as temp:
            case_dir = Path(temp) / "case"
            result = build_and_launch(
                electro_selectors={
                    "myocardiumSolver": "singleCellSolver",
                    "ionicModel": "AlievPanfilov",
                    "tissue": "myocyte",
                },
                physics_selectors={"type": "electroModel"},
                case_dir=case_dir,
            )
            self.assertEqual(result["status"], "complete")
            self.assertEqual(len(result["results"]), 1)
            self.assertEqual(result["results"][0]["status"], "ok")

            # The engine writes manifests under <case_dir>/output by
            # generic_case convention; locate run_manifest.
            output_dir = case_dir / "output"
            run_manifest_path = output_dir / "run_manifest.json"
            self.assertTrue(run_manifest_path.exists())

            run_manifest = json.loads(run_manifest_path.read_text())
            self.assertEqual(run_manifest["status"], "completed")
            self.assertIsNotNone(run_manifest["artifacts_realized_path"])

            realized = json.loads(
                Path(run_manifest["artifacts_realized_path"]).read_text()
            )
            self.assertEqual(realized["schema_version"], "1.1")
            self.assertGreaterEqual(len(realized["cases"]), 1)
            total_matched = sum(c["matched_count"] for c in realized["cases"])
            self.assertGreater(
                total_matched, 0,
                "expected at least one matched artifact after a real cardiacFoam run",
            )


if __name__ == "__main__":
    unittest.main()
