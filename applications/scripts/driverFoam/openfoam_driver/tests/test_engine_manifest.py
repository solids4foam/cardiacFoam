from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.engine import DriverEngine
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text())


class TestDriverEngineManifest(unittest.TestCase):
    def test_run_simulations_writes_rich_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            applied_cases: list[str] = []
            ran_cases: list[str] = []

            def build_cases() -> list[CaseConfig]:
                return [
                    CaseConfig("caseA", {"value": 1}),
                    CaseConfig("caseB", {"value": 2}),
                ]

            def apply_case(_: Path, case: CaseConfig) -> None:
                applied_cases.append(case.case_id)

            def run_case(_: Path, __: Path, case: CaseConfig) -> None:
                ran_cases.append(case.case_id)

            spec = TutorialSpec(
                name="dummy",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=build_cases,
                apply_case=apply_case,
                run_case=run_case,
            )

            engine = DriverEngine(spec=spec, requested_action="all")
            results = engine.run_simulations()

            manifest = _load_json(output_dir / "run_manifest.json")
            report_path = output_dir / "run_report.md"
            self.assertEqual(applied_cases, ["caseA", "caseB"])
            self.assertEqual(ran_cases, ["caseA", "caseB"])
            self.assertEqual(len(results), 2)
            self.assertEqual(manifest["schema_version"], "2.0")
            self.assertEqual(manifest["requested_action"], "all")
            self.assertEqual(manifest["status"], "completed")
            self.assertEqual(manifest["postprocess_status"], "not_started")
            self.assertEqual(manifest["completed_cases"], 2)
            self.assertEqual(manifest["failed_cases"], 0)
            self.assertIsNotNone(manifest["run_id"])
            self.assertIsNotNone(manifest["started_at_utc"])
            self.assertIsNotNone(manifest["finished_at_utc"])
            self.assertEqual(manifest["human_report_path"], str(report_path))
            self.assertEqual(manifest["results"][0]["index"], 1)
            self.assertEqual(manifest["results"][1]["total_cases"], 2)
            self.assertTrue(report_path.exists())
            report_text = report_path.read_text()
            self.assertIn("# driverFoam Run Report", report_text)
            self.assertIn("### caseA", report_text)
            self.assertIn("### caseB", report_text)

    def test_run_postprocess_writes_manifest_and_plots_path(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            def postprocess(_: Path, out_dir: Path) -> None:
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / "plots.json").write_text('{"plots": []}')

            spec = TutorialSpec(
                name="dummy",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=lambda: [],
                apply_case=lambda _case_root, _case: None,
                run_case=lambda _case_root, _setup_root, _case: None,
                postprocess=postprocess,
            )

            engine = DriverEngine(spec=spec, requested_action="post")
            engine.run_postprocess()

            manifest = _load_json(output_dir / "run_manifest.json")
            self.assertEqual(manifest["requested_action"], "post")
            self.assertEqual(manifest["status"], "completed")
            self.assertEqual(manifest["postprocess_status"], "completed")
            self.assertEqual(manifest["total_cases"], 0)
            self.assertEqual(manifest["results"], [])
            self.assertEqual(
                manifest["plots_manifest_path"],
                str(output_dir / "plots.json"),
            )
            self.assertTrue((output_dir / "run_report.md").exists())

    def test_run_postprocess_collects_outputs_before_postprocessing(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            collected: list[str] = []
            seen_in_postprocess: list[str] = []

            def collect_outputs(case_dir: Path, out_dir: Path) -> None:
                self.assertEqual(case_dir, case_root)
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / "3D_80_cells_implicit.dat").write_text("fresh-output")
                collected.append("ok")

            def postprocess(_: Path, out_dir: Path) -> None:
                seen_in_postprocess.extend(path.name for path in out_dir.glob("*.dat"))
                (out_dir / "plots.json").write_text('{"plots": []}')

            spec = TutorialSpec(
                name="dummy",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=lambda: [],
                apply_case=lambda _case_root, _case: None,
                run_case=lambda _case_root, _setup_root, _case: None,
                collect_outputs=collect_outputs,
                postprocess=postprocess,
            )

            engine = DriverEngine(spec=spec, requested_action="post")
            engine.run_postprocess()

            self.assertEqual(collected, ["ok"])
            self.assertEqual(seen_in_postprocess, ["3D_80_cells_implicit.dat"])
            self.assertTrue((output_dir / "3D_80_cells_implicit.dat").exists())


if __name__ == "__main__":
    unittest.main()
