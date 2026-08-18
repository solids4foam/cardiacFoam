"""The strict run payload must carry an artifact reconciliation report."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


def _artifact(**kwargs):
    from openfoam_driver.core.runtime.models import DataArtifact
    defaults = {"artifact_id": "x", "path_pattern": "out/x.dat", "format": "csv_probe"}
    defaults.update(kwargs)
    return DataArtifact(**defaults)


class TestReconciliationPayload(unittest.TestCase):
    def test_reports_matched_artifact(self) -> None:
        from openfoam_driver.cli import _reconciliation_payload
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp)
            (case_root / "out").mkdir()
            (case_root / "out" / "x.dat").write_text("data")
            payload = _reconciliation_payload(case_root, (_artifact(),))
            self.assertEqual(payload["matched_count"], 1)
            self.assertEqual(payload["missing_count"], 0)

    def test_reports_missing_artifact(self) -> None:
        from openfoam_driver.cli import _reconciliation_payload
        with tempfile.TemporaryDirectory() as tmp:
            payload = _reconciliation_payload(Path(tmp), (_artifact(),))
            self.assertEqual(payload["matched_count"], 0)
            self.assertEqual(payload["missing_count"], 1)

    def test_empty_expected_artifacts_is_safe(self) -> None:
        from openfoam_driver.cli import _reconciliation_payload
        with tempfile.TemporaryDirectory() as tmp:
            payload = _reconciliation_payload(Path(tmp), None)
            self.assertEqual(payload["predicted_count"], 0)


if __name__ == "__main__":
    unittest.main()
