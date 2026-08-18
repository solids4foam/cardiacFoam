"""Tolerances are transcribed from committed references, never chosen post hoc."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


REFERENCE_TEXT = """\
# file time variable expected tolerance
postProcessing/probes/0/Vm.dat 1.5 P1 -80.1 0.5
postProcessing/probes/0/Vm.dat 2.0 P2 -79.4 0.5
"""


class TestTranscription(unittest.TestCase):
    def test_one_row_per_reference_point(self) -> None:
        from openfoam_driver.tests.equivalence.protocol import transcribe_reference
        rows = transcribe_reference(
            "electrophysiologyProtocols/singleCell",
            "regression/singleCell.reference",
            REFERENCE_TEXT,
        )
        self.assertEqual(len(rows), 2)

    def test_tolerance_and_provenance_are_carried(self) -> None:
        from openfoam_driver.tests.equivalence.protocol import transcribe_reference
        rows = transcribe_reference(
            "electrophysiologyProtocols/singleCell",
            "regression/singleCell.reference",
            REFERENCE_TEXT,
        )
        row = rows[0]
        self.assertEqual(row.tolerance, 0.5)
        self.assertEqual(row.variable, "P1")
        self.assertEqual(row.source_reference, "regression/singleCell.reference")
        self.assertIn("committed", row.rationale)

    def test_non_columnar_reference_yields_no_rows(self) -> None:
        from openfoam_driver.tests.equivalence.protocol import transcribe_reference
        rows = transcribe_reference("x", "y.reference", "kind key metric a b c d\n")
        self.assertEqual(rows, [])


METRIC_TEXT = """\
# kind       key        metric      expected          tolerance
summary     cells      value       20                0
error       Vm         L1          0.000293945       1e-7
error       Vm         Linf        0.00148244        1e-7
topology    pvjNodes   count       8                 0
"""


class TestMetricTranscription(unittest.TestCase):
    """The MMS references carry the solver-emitted norm tolerances.

    They are already committed, so they need no pilot -- only a parser for the
    `kind key metric expected tolerance` layout the columnar parser rejects.
    """

    def test_transcribes_every_metric_row(self) -> None:
        from openfoam_driver.tests.equivalence.protocol import (
            transcribe_metric_reference,
        )
        rows = transcribe_metric_reference(
            "manufacturedSolutions/bidomain", "r.reference", METRIC_TEXT
        )
        self.assertEqual(len(rows), 4)

    def test_captures_norm_tolerances(self) -> None:
        from openfoam_driver.tests.equivalence.protocol import (
            transcribe_metric_reference,
        )
        rows = transcribe_metric_reference("c", "r.reference", METRIC_TEXT)
        norms = {(r.key, r.metric): r.tolerance for r in rows if r.kind == "error"}
        self.assertEqual(norms[("Vm", "L1")], 1e-7)
        self.assertEqual(norms[("Vm", "Linf")], 1e-7)

    def test_columnar_reference_yields_no_metric_rows(self) -> None:
        from openfoam_driver.tests.equivalence.protocol import (
            transcribe_metric_reference,
        )
        self.assertEqual(transcribe_metric_reference("c", "r", REFERENCE_TEXT), [])

    def test_metric_reference_yields_no_columnar_rows(self) -> None:
        from openfoam_driver.tests.equivalence.protocol import transcribe_reference
        self.assertEqual(transcribe_reference("c", "r", METRIC_TEXT), [])


class TestRoundTrip(unittest.TestCase):
    def test_written_protocol_loads_back_identically(self) -> None:
        from openfoam_driver.tests.equivalence.protocol import (
            transcribe_reference, transcribe_metric_reference,
            write_protocol, load_protocol,
        )
        rows = transcribe_reference("case/a", "r.reference", REFERENCE_TEXT)
        metric_rows = transcribe_metric_reference("case/b", "m.reference", METRIC_TEXT)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "equivalence_protocol.yaml"
            write_protocol(rows, path, metric_rows=metric_rows)
            loaded = load_protocol(path)
            # write_protocol sorts deliberately, so compare content not order.
            self.assertEqual(set(loaded.rows), set(rows))
            self.assertEqual(set(loaded.metric_rows), set(metric_rows))


if __name__ == "__main__":
    unittest.main()
