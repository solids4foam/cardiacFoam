"""Contract tests for the DataArtifact vocabulary (plan v2 phase 2).

DataArtifact is the shared output-description language between the engine
(run-side predictions written to artifacts_manifest.json) and the utility
catalog (static declarations in utility.manifest.toml). The fields, defaults,
and ArtifactFormat enum below are part of the agent-facing contract — every
change here is observed by downstream consumers.
"""
from __future__ import annotations

import dataclasses
import typing
import unittest

from openfoam_driver.core.runtime.models import ArtifactFormat, DataArtifact


class TestDataArtifact(unittest.TestCase):
    def test_constructs_with_required_fields_only(self) -> None:
        artifact = DataArtifact(
            artifact_id="vm_probe",
            path_pattern="postProcessing/probes/{time}/Vm",
            format="csv_probe",
        )
        self.assertEqual(artifact.artifact_id, "vm_probe")
        self.assertEqual(artifact.path_pattern, "postProcessing/probes/{time}/Vm")
        self.assertEqual(artifact.format, "csv_probe")

    def test_defaults_are_safe_for_predictor_merging(self) -> None:
        """Defaults must let predict_data_artifacts merge static + derived
        artifacts without None-vs-tuple ambiguity (plan section 2.1)."""
        artifact = DataArtifact(
            artifact_id="x",
            path_pattern="foo",
            format="openfoam_log",
        )
        self.assertEqual(artifact.variables, ())  # never None
        self.assertEqual(artifact.description, "")
        self.assertEqual(artifact.produced_by, "")
        self.assertIs(artifact.optional, False)
        self.assertIs(artifact.time_indexed, False)

    def test_is_frozen(self) -> None:
        """Artifacts are value objects embedded in agent manifests; mutation
        would silently desync the manifest from later reads."""
        artifact = DataArtifact(
            artifact_id="x",
            path_pattern="foo",
            format="openfoam_log",
        )
        with self.assertRaises(dataclasses.FrozenInstanceError):
            artifact.artifact_id = "y"  # type: ignore[misc]

    def test_accepts_variables_tuple(self) -> None:
        artifact = DataArtifact(
            artifact_id="ionic",
            path_pattern="postProcessing/cellModel.dat",
            format="csv_sweep",
            variables=("Vm", "Iion", "Cai"),
        )
        self.assertEqual(artifact.variables, ("Vm", "Iion", "Cai"))


class TestArtifactFormatLiteral(unittest.TestCase):
    """ArtifactFormat is a closed Literal — adding values requires updating
    every consumer that branches on format. Lock the set here."""

    _EXPECTED_FORMATS: frozenset[str] = frozenset({
        "csv_probe",
        "csv_sweep",
        "vtk_sequence",
        "openfoam_time_dirs",
        "openfoam_log",
        "json_summary",
    })

    def test_literal_values_match_documented_set(self) -> None:
        actual = frozenset(typing.get_args(ArtifactFormat))
        self.assertEqual(
            actual, self._EXPECTED_FORMATS,
            "ArtifactFormat enum changed — update plan section 2.1 and every "
            "consumer that branches on format before changing this assertion.",
        )


if __name__ == "__main__":
    unittest.main()
