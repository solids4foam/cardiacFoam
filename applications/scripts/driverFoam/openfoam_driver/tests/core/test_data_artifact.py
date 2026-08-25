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
#     test_data_artifact
#
# Description
#     Tests data artifact logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Contract tests for the DataArtifact vocabulary (plan v2 phase 2).

DataArtifact is the shared output-description language between the engine
and the utility
catalog (static declarations in utility.manifest.toml). The fields, defaults,
and ArtifactFormat enum below are part of the agent-facing contract — every
change here is observed by downstream consumers.
"""
from __future__ import annotations

import dataclasses
import typing
import unittest

from openfoam_driver.core.runtime.models import (
    ArtifactFormat,
    DataArtifact,
    expand_path_pattern,
)


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

    def test_construction_rejects_unknown_placeholder(self) -> None:
        """A typo (e.g. {caseId}) in path_pattern must fail at construction,
        not silently propagate where it can
        only be detected when an agent tries to expand it later."""
        with self.assertRaises(ValueError) as ctx:
            DataArtifact(
                artifact_id="typo",
                path_pattern="postProcessing/{caseId}.txt",
                format="csv_probe",
            )
        self.assertIn("caseId", str(ctx.exception))

    def test_construction_with_known_placeholders_succeeds(self) -> None:
        DataArtifact(
            artifact_id="ok",
            path_pattern="results/{case_id}/{time}/Vm",
            format="openfoam_time_dirs",
        )

    def test_construction_with_no_placeholders_succeeds(self) -> None:
        DataArtifact(
            artifact_id="static",
            path_pattern="postProcessing/exact_error.json",
            format="json_summary",
        )

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


class TestExpandPathPattern(unittest.TestCase):
    """Path-pattern placeholders ({case_id}, {time}) are the only documented
    substitution language. The helper enforces the closed set so a predictor
    or utility-manifest author cannot silently invent a third placeholder."""

    def test_no_placeholders_returns_pattern_unchanged(self) -> None:
        out = expand_path_pattern("postProcessing/probes.dat", case_id="c1")
        self.assertEqual(out, "postProcessing/probes.dat")

    def test_substitutes_case_id(self) -> None:
        out = expand_path_pattern("results/{case_id}/log", case_id="caseA")
        self.assertEqual(out, "results/caseA/log")

    def test_substitutes_time(self) -> None:
        out = expand_path_pattern(
            "postProcessing/probes/{time}/Vm", case_id="c1", time="0.01"
        )
        self.assertEqual(out, "postProcessing/probes/0.01/Vm")

    def test_substitutes_both_placeholders(self) -> None:
        out = expand_path_pattern(
            "out/{case_id}/{time}/field", case_id="c2", time="0.5"
        )
        self.assertEqual(out, "out/c2/0.5/field")

    def test_repeated_placeholder_is_substituted_everywhere(self) -> None:
        out = expand_path_pattern("{case_id}/{case_id}.log", case_id="x")
        self.assertEqual(out, "x/x.log")

    def test_missing_case_id_for_pattern_that_needs_it_raises(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            expand_path_pattern("results/{case_id}/log", case_id=None)
        self.assertIn("case_id", str(ctx.exception))

    def test_missing_time_for_pattern_that_needs_it_raises(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            expand_path_pattern(
                "postProcessing/{time}/Vm", case_id="c1", time=None
            )
        self.assertIn("time", str(ctx.exception))

    def test_unknown_placeholder_raises(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            expand_path_pattern("out/{foo}/x", case_id="c1")
        self.assertIn("foo", str(ctx.exception))
        self.assertIn("unknown placeholder", str(ctx.exception).lower())

    def test_unused_kwargs_are_tolerated(self) -> None:
        """Passing time= when the pattern has no {time} is not an error.
        Callers compose artifacts uniformly; they should not have to inspect
        each pattern before calling."""
        out = expand_path_pattern(
            "postProcessing/static.csv", case_id="c1", time="0.01"
        )
        self.assertEqual(out, "postProcessing/static.csv")


if __name__ == "__main__":
    unittest.main()
