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
#     test_table_writer
#
# Description
#     Tests table writer logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import tempfile
import unittest
from datetime import datetime
from pathlib import Path

from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter


class TestTableWriter(unittest.TestCase):
    def test_writes_csv_with_envelope(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            rows = [
                {"case_id": "case_A", "DX_mm": 0.1, "activation_ms": 42.3},
                {"case_id": "case_B", "DX_mm": 0.2, "activation_ms": 55.1},
            ]
            meta = TableMetadata(
                entry="TestTutorial",
                units={"activation_ms": "ms", "DX_mm": "mm"},
            )
            artifacts = TableWriter.write(rows, output_dir, "test_summary", "Test label", meta)

            csv_path = output_dir / "test_summary.csv"
            self.assertTrue(csv_path.exists())
            text = csv_path.read_text()
            self.assertIn("# entry: TestTutorial", text)
            self.assertIn("# generated_at:", text)
            self.assertIn('"activation_ms": "ms"', text)
            self.assertIn("case_id,DX_mm,activation_ms", text)
            self.assertIn("case_A,0.1,42.3", text)
            self.assertIn("case_B,0.2,55.1", text)

    def test_writes_html_with_metadata_and_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            rows = [{"col_a": "x", "col_b": 1}]
            meta = TableMetadata(entry="HtmlTest", units={"col_b": "ms"})
            TableWriter.write(rows, output_dir, "html_test", "HTML label", meta)

            html_text = (output_dir / "html_test.html").read_text()
            self.assertIn("HtmlTest", html_text)
            self.assertIn("col_a", html_text)
            self.assertIn("<td>x</td>", html_text)

    def test_returns_two_artifacts_with_correct_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            meta = TableMetadata(entry="ArtifactTest", units={})
            artifacts = TableWriter.write(
                [{"x": 1}], output_dir, "art_stem", "Art label", meta
            )
            self.assertEqual(len(artifacts), 2)
            by_format = {a["format"]: a for a in artifacts}
            self.assertIn("csv", by_format)
            self.assertIn("html", by_format)
            for a in artifacts:
                self.assertEqual(a["kind"], "table")
                self.assertEqual(a["label"], "Art label")
                self.assertIn("path", a)

    def test_handles_empty_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            meta = TableMetadata(entry="EmptyTest", units={})
            artifacts = TableWriter.write([], output_dir, "empty_stem", "Empty", meta)
            self.assertEqual(len(artifacts), 2)
            csv_text = (output_dir / "empty_stem.csv").read_text()
            self.assertIn("# entry: EmptyTest", csv_text)

    def test_metadata_autofills_generated_at(self) -> None:
        meta = TableMetadata(entry="T", units={})
        self.assertNotEqual(meta.generated_at, "")
        # Must be a parseable ISO-8601 string
        datetime.fromisoformat(meta.generated_at)

    def test_metadata_preserves_explicit_generated_at(self) -> None:
        meta = TableMetadata(entry="T", units={}, generated_at="2026-01-01T00:00:00+00:00")
        self.assertEqual(meta.generated_at, "2026-01-01T00:00:00+00:00")

    def test_artifact_paths_are_relative_filenames(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            meta = TableMetadata(entry="T", units={})
            artifacts = TableWriter.write([{"a": 1}], output_dir, "stem", "L", meta)
            for a in artifacts:
                # path must be just the filename, not absolute
                self.assertFalse(Path(a["path"]).is_absolute())
                self.assertIn("stem", a["path"])


if __name__ == "__main__":
    unittest.main()
