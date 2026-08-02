"""Tests for data_artifact_from_json."""
from __future__ import annotations

import unittest

from openfoam_driver.core.runtime.models import DataArtifact, data_artifact_from_json


class TestDataArtifactFromJson(unittest.TestCase):
    def test_round_trips_a_full_artifact(self) -> None:
        original = DataArtifact(
            artifact_id="vm_probe",
            path_pattern="postProcessing/{case_id}/Vm.csv",
            format="csv_probe",
            variables=("t", "Vm"),
            description="membrane potential probe",
            produced_by="cardiacFoam",
            optional=True,
            time_indexed=False,
        )
        from dataclasses import asdict
        payload = asdict(original)
        payload["variables"] = list(original.variables)
        rebuilt = data_artifact_from_json(payload)
        self.assertEqual(rebuilt, original)

    def test_minimal_payload_uses_defaults(self) -> None:
        rebuilt = data_artifact_from_json({
            "artifact_id": "x",
            "path_pattern": "a/b.dat",
            "format": "openfoam_log",
        })
        self.assertEqual(rebuilt.variables, ())
        self.assertEqual(rebuilt.description, "")
        self.assertEqual(rebuilt.produced_by, "")
        self.assertFalse(rebuilt.optional)
        self.assertFalse(rebuilt.time_indexed)

    def test_variables_is_coerced_to_tuple_of_str(self) -> None:
        rebuilt = data_artifact_from_json({
            "artifact_id": "x",
            "path_pattern": "a.dat",
            "format": "csv_probe",
            "variables": ["t", "Vm"],
        })
        self.assertIsInstance(rebuilt.variables, tuple)
        self.assertEqual(rebuilt.variables, ("t", "Vm"))

    def test_unknown_placeholder_still_raises_via_post_init(self) -> None:
        with self.assertRaises(ValueError):
            data_artifact_from_json({
                "artifact_id": "x",
                "path_pattern": "a/{bogus}.dat",
                "format": "csv_probe",
            })


if __name__ == "__main__":
    unittest.main()
