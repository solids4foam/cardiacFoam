from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.postprocessing.single_cell_mode_compare import (
    compare_single_cell_modes,
    _mode_sort_key,
)


def _write_manifest(
    root: Path,
    mode: str,
    ionic_model: str,
    *,
    status: str = "ok",
) -> Path:
    output_dir = root / mode / "postProcessing"
    output_dir.mkdir(parents=True)
    manifest = output_dir / "run_manifest.json"
    manifest.write_text(
        "{\n"
        '  "schema_version": "2.3",\n'
        '  "output_dir": "' + str(output_dir) + '",\n'
        '  "results": [\n'
        "    {\n"
        '      "case_id": "' + ionic_model + '_myocyte",\n'
        '      "status": "' + status + '",\n'
        '      "duration_s": 0.0,\n'
        '      "params": {"ionicModel": "' + ionic_model + '", "tissue": "myocyte"}\n'
        "    }\n"
        "  ]\n"
        "}\n"
    )
    return manifest


def _write_trace(path: Path, rows: list[tuple[float, float, float]]) -> None:
    path.write_text(
        "time Vm gate\n"
        + "".join(f"{time} {vm} {gate}\n" for time, vm, gate in rows)
    )


class TestSingleCellModeCompare(unittest.TestCase):
    def test_compares_batched_trace_against_cpu_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cpu_manifest = _write_manifest(root, "cpu", "Gaur")
            batched_manifest = _write_manifest(root, "batched_rl", "GaurcompactBatched")

            _write_trace(
                cpu_manifest.parent / "Gaur_myocyte_S1_1000.txt",
                [(0.0, -80.0, 0.1), (0.5, -70.0, 0.2), (1.0, -60.0, 0.3)],
            )
            _write_trace(
                batched_manifest.parent / "GaurcompactBatched_myocyte_S1_1000.txt",
                [(0.0, -79.0, 0.1), (1.0, -58.0, 0.4)],
            )

            output_dir = root / "comparison"
            result = compare_single_cell_modes(
                {"cpu": cpu_manifest, "batched_rl": batched_manifest},
                output_dir,
                make_plots=False,
            )

            self.assertTrue(result["artifacts"])
            summary = output_dir / "single_cell_average_differences.csv"
            self.assertTrue(summary.exists())

            with summary.open() as handle:
                rows = list(csv.DictReader(handle))

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["model"], "Gaur")
            self.assertEqual(rows[0]["tissue"], "myocyte")
            self.assertEqual(rows[0]["mode"], "batched_rl")
            self.assertEqual(rows[0]["variable"], "Vm")
            self.assertEqual(rows[0]["status"], "ok")
            self.assertAlmostEqual(float(rows[0]["mean_abs_diff"]), 1.5)

    def test_skips_failed_manifest_entry_even_when_trace_exists(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cpu_manifest = _write_manifest(root, "cpu", "Gaur")
            failed_manifest = _write_manifest(
                root,
                "rl",
                "GaurcompactBatched",
                status="failed",
            )

            _write_trace(
                cpu_manifest.parent / "Gaur_myocyte_S1_1000.txt",
                [(0.0, -80.0, 0.1), (1.0, -60.0, 0.3)],
            )
            _write_trace(
                failed_manifest.parent / "GaurcompactBatched_myocyte_S1_1000.txt",
                [(0.0, -79.0, 0.1), (1.0, -58.0, 0.4)],
            )

            output_dir = root / "comparison"
            result = compare_single_cell_modes(
                {"cpu": cpu_manifest, "rl": failed_manifest},
                output_dir,
                make_plots=False,
            )

            self.assertIn(
                "rl: skipping case GaurcompactBatched_myocyte with status failed",
                result["warnings"],
            )

            with (output_dir / "single_cell_average_differences.csv").open() as handle:
                rows = list(csv.DictReader(handle))

            self.assertEqual(rows, [])

    def test_step_modes_sort_by_family_then_descending_step_count(self) -> None:
        modes = [
            "soa_step010",
            "euler_step020",
            "rl_step030",
            "cpu",
            "euler_step050",
            "soa_step050",
            "rl_step010",
            "euler_step040",
        ]

        self.assertEqual(
            sorted(modes, key=_mode_sort_key),
            [
                "cpu",
                "euler_step050",
                "euler_step040",
                "euler_step020",
                "rl_step030",
                "rl_step010",
                "soa_step050",
                "soa_step010",
            ],
        )


if __name__ == "__main__":
    unittest.main()
