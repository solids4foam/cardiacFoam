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
#     test_reconciler
#
# Description
#     Tests reconciler logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the post-run artifact reconciler.

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


class TestRealPredictorOutputShapes(unittest.TestCase):
    """Regression tests for the 2026-05-21 audit findings:

    - PDE solver predictors emit `path_pattern="{time}"` — the OpenFOAM
      time directory itself, not a file inside it. The reconciler must
      accept directories as matches.
    - The single-cell predictor emits `path_pattern="postProcessing/{case_id}.txt"`.
      With no case_id supplied, `{case_id}` becomes a glob wildcard so the
      reconciler still finds per-case outputs across a sweep.
    """

    def test_time_directory_alone_is_a_match(self) -> None:
        """`path_pattern="{time}"` (the monodomain predictor's actual
        output) must match the time directory itself, reported as kind='dir'."""
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            for t in ("0", "0.001", "0.002"):
                (case_root / t).mkdir()
                # Each time dir has a Vm file inside — typical OpenFOAM.
                (case_root / t / "Vm").write_bytes(b"v")
            (case_root / "constant").mkdir()  # must be ignored

            artifact = _make_artifact(
                artifact_id="myocardium_time_series",
                path_pattern="{time}",
                format="openfoam_time_dirs",
                time_indexed=True,
            )
            report = reconcile_artifacts(case_root, (artifact,))
            entry = report.artifacts[0]
            self.assertEqual(entry["status"], "matched")
            self.assertEqual(len(entry["matched_files"]), 3)
            for match in entry["matched_files"]:
                self.assertEqual(match["kind"], "dir")
                self.assertIn("entries", match)
                self.assertGreaterEqual(match["entries"], 1)

    def test_case_id_glob_matches_every_per_case_file(self) -> None:
        """The single-cell predictor's `postProcessing/{case_id}.txt`
        with no case_id supplied must glob across every existing per-case
        file in postProcessing/."""
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            (case_root / "postProcessing").mkdir()
            for name in ("TNNP_M.txt", "TNNP_epi.txt", "AlievPanfilov_myocyte.txt"):
                (case_root / "postProcessing" / name).write_bytes(b"x")

            artifact = _make_artifact(
                artifact_id="single_cell_trace",
                path_pattern="postProcessing/{case_id}.txt",
                format="csv_sweep",
            )
            report = reconcile_artifacts(case_root, (artifact,))
            entry = report.artifacts[0]
            self.assertEqual(entry["status"], "matched")
            self.assertEqual(len(entry["matched_files"]), 3)
            for match in entry["matched_files"]:
                self.assertEqual(match["kind"], "file")

    def test_case_id_literal_substitution_when_supplied(self) -> None:
        """When `case_id="TNNP_M"` is supplied, only the matching per-case
        file is matched — not the others."""
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            (case_root / "postProcessing").mkdir()
            for name in ("TNNP_M.txt", "TNNP_epi.txt", "AlievPanfilov_myocyte.txt"):
                (case_root / "postProcessing" / name).write_bytes(b"x")

            artifact = _make_artifact(
                artifact_id="single_cell_trace",
                path_pattern="postProcessing/{case_id}.txt",
                format="csv_sweep",
            )
            report = reconcile_artifacts(
                case_root, (artifact,), case_id="TNNP_M",
            )
            entry = report.artifacts[0]
            self.assertEqual(entry["status"], "matched")
            self.assertEqual(len(entry["matched_files"]), 1)
            self.assertTrue(entry["matched_files"][0]["path"].endswith("TNNP_M.txt"))

    def test_case_id_literal_missing_when_no_matching_file(self) -> None:
        """case_id substituted literally but no matching file → missing."""
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            (case_root / "postProcessing").mkdir()
            (case_root / "postProcessing" / "TNNP_M.txt").write_bytes(b"x")

            artifact = _make_artifact(
                artifact_id="single_cell_trace",
                path_pattern="postProcessing/{case_id}.txt",
                format="csv_sweep",
            )
            report = reconcile_artifacts(
                case_root, (artifact,), case_id="not_a_real_case",
            )
            self.assertEqual(report.artifacts[0]["status"], "missing")

    def test_case_id_field_on_report(self) -> None:
        """The report carries the case_id used (or None for whole-run)."""
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            artifact = _make_artifact(path_pattern="absent.dat")
            report1 = reconcile_artifacts(case_root, (artifact,))
            report2 = reconcile_artifacts(case_root, (artifact,), case_id="c1")
            self.assertIsNone(report1.case_id)
            self.assertEqual(report2.case_id, "c1")


class TestReconcilerCaseId(unittest.TestCase):
    """The reconciler accepts an optional case_id label so per-case
    reports can be attributed in multi-case sweeps."""

    def test_case_id_is_recorded_on_report(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            artifact = _make_artifact(
                artifact_id="x", path_pattern="absent.dat", format="csv_probe",
            )
            report = reconcile_artifacts(case_root, (artifact,), case_id="myCase")
            self.assertEqual(report.case_id, "myCase")

    def test_default_case_id_is_none(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp)
            artifact = _make_artifact(
                artifact_id="x", path_pattern="absent.dat", format="csv_probe",
            )
            report = reconcile_artifacts(case_root, (artifact,))
            self.assertIsNone(report.case_id)


class TestFileHashing(unittest.TestCase):
    def test_matched_file_entry_carries_sha256(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp)
            (case_root / "out").mkdir()
            (case_root / "out" / "x.dat").write_text("hello")
            report = reconcile_artifacts(case_root, (_make_artifact(),))
            entry = report.artifacts[0]
            self.assertEqual(entry["status"], "matched")
            match = entry["matched_files"][0]
            # sha256("hello")
            self.assertEqual(
                match["sha256"],
                "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824",
            )

    def test_directory_entry_has_no_sha256(self) -> None:
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp)
            (case_root / "out").mkdir()
            (case_root / "out" / "d").mkdir()
            report = reconcile_artifacts(
                case_root, (_make_artifact(path_pattern="out/d"),)
            )
            match = report.artifacts[0]["matched_files"][0]
            self.assertEqual(match["kind"], "dir")
            self.assertNotIn("sha256", match)


class TestReportSerialization(unittest.TestCase):
    def test_to_json_round_trips_summary_fields(self) -> None:
        import json
        from openfoam_driver.core.runtime.reconciler import reconcile_artifacts
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp)
            report = reconcile_artifacts(
                case_root, (_make_artifact(),), case_id="c1"
            )
            payload = report.to_json()
            self.assertEqual(payload["case_id"], "c1")
            self.assertEqual(payload["predicted_count"], 1)
            self.assertEqual(payload["missing_count"], 1)
            self.assertIsInstance(payload["artifacts"], list)
            json.dumps(payload)  # must be JSON-serializable


if __name__ == "__main__":
    unittest.main()
