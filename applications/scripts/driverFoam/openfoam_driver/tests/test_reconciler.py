"""Tests for the post-run artifact reconciler (plan §10).

The reconciler walks a list of predicted `DataArtifact` objects and
compares them against the on-disk state of `case_root`. Output is a
report that classifies each artifact as matched / missing, and lists
the actual files that matched (for matched artifacts).
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


def _make_artifact(**kwargs):
    """Tiny builder so tests can spell out only the fields they care about."""
    from openfoam_driver.core.runtime.models import DataArtifact
    defaults = {
        "artifact_id": "x",
        "path_pattern": "out/x.dat",
        "format": "csv_probe",
    }
    defaults.update(kwargs)
    return DataArtifact(**defaults)


class TestReconcilerModule(unittest.TestCase):
    def test_module_exports_reconcile_artifacts(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        self.assertTrue(callable(reconcile_artifacts))

    def test_module_exports_report_dataclass(self) -> None:
        from openfoam_driver.core.runtime.reconciler import ReconciliationReport
        self.assertTrue(hasattr(ReconciliationReport, "__dataclass_fields__"))


class TestNonTimeIndexedReconciliation(unittest.TestCase):
    """A non-time-indexed artifact has a single expected path; reconciler
    checks the file exists and reports its size."""

    def test_present_file_is_matched(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            (case_root / "postProcessing").mkdir()
            file_path = case_root / "postProcessing" / "probes.dat"
            file_path.write_bytes(b"abcdef")

            artifact = _make_artifact(
                artifact_id="probes",
                path_pattern="postProcessing/probes.dat",
                format="csv_probe",
            )
            report = reconcile_artifacts(case_root, (artifact,))
            self.assertEqual(len(report.artifacts), 1)
            entry = report.artifacts[0]
            self.assertEqual(entry["artifact_id"], "probes")
            self.assertEqual(entry["status"], "matched")
            self.assertEqual(len(entry["matched_files"]), 1)
            self.assertEqual(entry["matched_files"][0]["size_bytes"], 6)
            self.assertTrue(entry["matched_files"][0]["path"].endswith("probes.dat"))

    def test_missing_file_is_classified_missing(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            artifact = _make_artifact(
                artifact_id="absent",
                path_pattern="postProcessing/absent.dat",
                format="csv_probe",
            )
            report = reconcile_artifacts(case_root, (artifact,))
            entry = report.artifacts[0]
            self.assertEqual(entry["status"], "missing")
            self.assertEqual(entry["matched_files"], [])

    def test_missing_optional_artifact_is_flagged_as_optional(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            artifact = _make_artifact(
                artifact_id="opt",
                path_pattern="postProcessing/opt.dat",
                format="csv_probe",
                optional=True,
            )
            report = reconcile_artifacts(case_root, (artifact,))
            entry = report.artifacts[0]
            self.assertEqual(entry["status"], "missing")
            self.assertTrue(entry["optional"])


class TestTimeIndexedReconciliation(unittest.TestCase):
    """Time-indexed artifacts (e.g. OpenFOAM time directories) match
    against multiple on-disk files — one per time dir that exists."""

    def test_glob_collects_every_time_directory(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            for t in ("0", "0.001", "0.002"):
                (case_root / t).mkdir()
                (case_root / t / "Vm").write_bytes(b"v")
            (case_root / "constant").mkdir()
            (case_root / "constant" / "Vm").write_bytes(b"static")

            artifact = _make_artifact(
                artifact_id="vm_series",
                path_pattern="{time}/Vm",
                format="openfoam_time_dirs",
                time_indexed=True,
            )
            report = reconcile_artifacts(case_root, (artifact,))
            entry = report.artifacts[0]
            self.assertEqual(entry["status"], "matched")
            self.assertEqual(len(entry["matched_files"]), 3)
            matched_paths = sorted(m["path"] for m in entry["matched_files"])
            for expected_t in ("0/Vm", "0.001/Vm", "0.002/Vm"):
                self.assertTrue(
                    any(p.endswith(expected_t) for p in matched_paths),
                    f"missing match for {expected_t} in {matched_paths}",
                )

    def test_time_indexed_with_no_time_dirs_is_missing(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            artifact = _make_artifact(
                artifact_id="vm_series",
                path_pattern="{time}/Vm",
                format="openfoam_time_dirs",
                time_indexed=True,
            )
            report = reconcile_artifacts(case_root, (artifact,))
            self.assertEqual(report.artifacts[0]["status"], "missing")


class TestReportSummary(unittest.TestCase):
    def test_report_carries_case_root_and_predicted_count(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            artifacts = (
                _make_artifact(artifact_id="a1", path_pattern="x.dat", format="csv_probe"),
                _make_artifact(artifact_id="a2", path_pattern="y.dat", format="csv_probe"),
            )
            report = reconcile_artifacts(case_root, artifacts)
            self.assertEqual(report.case_root, str(case_root))
            self.assertEqual(len(report.artifacts), 2)
            self.assertEqual(report.predicted_count, 2)
            self.assertEqual(report.matched_count, 0)
            self.assertEqual(report.missing_count, 2)


if __name__ == "__main__":
    unittest.main()
