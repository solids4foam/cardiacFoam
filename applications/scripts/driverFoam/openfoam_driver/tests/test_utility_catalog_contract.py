"""Contract test: every utility directory has a manifest and the catalog is self-consistent.

When ``test_every_utility_dir_has_a_manifest`` fails, create a
``utility.manifest.toml`` in the missing utility directory.

When any other test fails, either the manifest TOML or the catalog loader has
drifted; fix the manifest file or the category constant in ``utility_catalog.py``.
"""

from __future__ import annotations

import textwrap
import unittest
import warnings
from pathlib import Path

from openfoam_driver.utility_catalog import (
    ALLOWED_CATEGORIES,
    ALLOWED_ARGUMENT_KINDS,
    MANIFEST_FILENAME,
    UTILITY_CATALOG,
    UtilityFlag,
    UtilityManifest,
    ProducesEntry,
    load_utility_manifests,
    _parse_manifest,
)

REPO_ROOT = Path(__file__).resolve().parents[5]
UTILITIES_ROOT = REPO_ROOT / "applications" / "utilities"


def _utility_dirs() -> list[Path]:
    """Return sorted list of direct subdirectories under UTILITIES_ROOT."""
    return sorted(p for p in UTILITIES_ROOT.iterdir() if p.is_dir())


def _write_minimal_toml(path: Path, extra: str = "") -> None:
    """Write a minimal valid manifest TOML to *path* (parent dir must exist)."""
    path.write_text(
        textwrap.dedent(
            f"""\
            name = "{path.parent.name}"
            description = "test utility"
            category = "verification"
            {extra}
            """
        ),
        encoding="utf-8",
    )


class TestUtilityCatalogContract(unittest.TestCase):
    def test_every_utility_dir_has_a_manifest(self) -> None:
        """Every directory under applications/utilities/ must have a manifest file."""
        missing = [
            d.name
            for d in _utility_dirs()
            if not (d / MANIFEST_FILENAME).exists()
        ]
        self.assertEqual(
            missing,
            [],
            f"Utility directories missing {MANIFEST_FILENAME}: {missing}",
        )

    def test_manifest_name_matches_directory(self) -> None:
        """Each loaded manifest's 'name' field must equal its parent directory name."""
        mismatches = []
        for name, manifest in UTILITY_CATALOG.items():
            dir_name = manifest.source_path.parent.name
            if manifest.name != dir_name:
                mismatches.append(
                    f"{manifest.source_path}: name={manifest.name!r} "
                    f"vs dir={dir_name!r}"
                )
        self.assertEqual(
            mismatches,
            [],
            f"name/directory mismatches: {mismatches}",
        )

    def test_no_duplicate_manifest_names(self) -> None:
        """The catalog size must equal the number of utility directories that have manifests."""
        dirs_with_manifests = [
            d for d in _utility_dirs() if (d / MANIFEST_FILENAME).exists()
        ]
        self.assertEqual(
            len(UTILITY_CATALOG),
            len(dirs_with_manifests),
            f"Catalog has {len(UTILITY_CATALOG)} entries but "
            f"{len(dirs_with_manifests)} manifest files exist; "
            "check for duplicate 'name' values across manifests.",
        )

    def test_required_fields_present(self) -> None:
        """Every manifest must have non-empty name, description, and category."""
        violations = []
        for name, manifest in UTILITY_CATALOG.items():
            if not manifest.name:
                violations.append(f"{name}: 'name' is empty")
            if not manifest.description:
                violations.append(f"{name}: 'description' is empty")
            if not manifest.category:
                violations.append(f"{name}: 'category' is empty")
        self.assertEqual(
            violations,
            [],
            f"Required-field violations: {violations}",
        )

    def test_categories_in_allowed_set(self) -> None:
        """Every manifest's category must be one of ALLOWED_CATEGORIES."""
        bad = [
            f"{name}: {manifest.category!r}"
            for name, manifest in UTILITY_CATALOG.items()
            if manifest.category not in ALLOWED_CATEGORIES
        ]
        self.assertEqual(
            bad,
            [],
            f"Categories not in ALLOWED_CATEGORIES {sorted(ALLOWED_CATEGORIES)}: {bad}",
        )


class TestArgumentKindEnum(unittest.TestCase):
    """Phase 4: argument_kind closed-enum validation in [[flags]] entries."""

    def _make_toml_with_flag_kind(self, kind: str) -> str:
        return textwrap.dedent(
            f"""\
            [[flags]]
            name = "-foo"
            description = "a flag"
            takes_value = true
            argument_kind = "{kind}"
            """
        )

    def test_valid_argument_kinds_accepted(self) -> None:
        """All six documented argument_kind values must parse without error."""
        import tempfile, os
        for kind in ALLOWED_ARGUMENT_KINDS:
            with tempfile.TemporaryDirectory() as tmp:
                d = Path(tmp) / "myutil"
                d.mkdir()
                p = d / "utility.manifest.toml"
                _write_minimal_toml(p, self._make_toml_with_flag_kind(kind))
                manifest = _parse_manifest(p)
                self.assertEqual(len(manifest.flags), 1)
                self.assertEqual(manifest.flags[0].argument_kind, kind)

    def test_invalid_argument_kind_rejected(self) -> None:
        """An unrecognised argument_kind value must raise ValueError at parse time."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "badutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, self._make_toml_with_flag_kind("integer"))
            with self.assertRaises(ValueError, msg="'integer' should not be a valid argument_kind"):
                _parse_manifest(p)

    def test_allowed_argument_kinds_set_is_closed(self) -> None:
        """ALLOWED_ARGUMENT_KINDS must equal exactly the six documented values."""
        self.assertEqual(
            ALLOWED_ARGUMENT_KINDS,
            frozenset({"scalar", "label", "path", "word", "word_list", "switch"}),
        )


class TestFlagNewOptionalFields(unittest.TestCase):
    """Phase 4: required and default optional fields on [[flags]]."""

    def _parse_with_flag(self, flag_toml: str) -> UtilityManifest:
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, flag_toml)
            return _parse_manifest(p)

    def test_required_defaults_to_false(self) -> None:
        """A flag with no 'required' field must default to False."""
        manifest = self._parse_with_flag(
            '[[flags]]\nname = "-x"\ndescription = "x"\ntakes_value = true\nargument_kind = "scalar"\n'
        )
        self.assertFalse(manifest.flags[0].required)

    def test_required_true_accepted(self) -> None:
        """A flag with required = true must parse as True."""
        manifest = self._parse_with_flag(
            '[[flags]]\nname = "-x"\ndescription = "x"\ntakes_value = true\nargument_kind = "scalar"\nrequired = true\n'
        )
        self.assertTrue(manifest.flags[0].required)

    def test_default_string_accepted(self) -> None:
        """A flag with a string default must parse and be accessible."""
        manifest = self._parse_with_flag(
            '[[flags]]\nname = "-x"\ndescription = "x"\ntakes_value = true\nargument_kind = "scalar"\ndefault = "3e-4"\n'
        )
        self.assertEqual(manifest.flags[0].default, "3e-4")

    def test_default_absent_is_none(self) -> None:
        """A flag with no default must yield None."""
        manifest = self._parse_with_flag(
            '[[flags]]\nname = "-x"\ndescription = "x"\ntakes_value = true\nargument_kind = "scalar"\n'
        )
        self.assertIsNone(manifest.flags[0].default)


class TestPositionalArgs(unittest.TestCase):
    """Phase 4: positional_args field on the manifest."""

    def test_positional_args_list_accepted(self) -> None:
        """positional_args with name, argument_kind, description must parse correctly."""
        import tempfile
        # Single-line inline tables are TOML 1.0 compatible.
        extra = textwrap.dedent(
            """\
            positional_args = [
              {name = "vtk_file", argument_kind = "path", description = "Input VTK file"},
              {name = "nsteps", argument_kind = "scalar", description = "Step count"},
            ]
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            manifest = _parse_manifest(p)
        self.assertEqual(len(manifest.positional_args), 2)
        first = manifest.positional_args[0]
        self.assertEqual(first.name, "vtk_file")
        self.assertEqual(first.argument_kind, "path")
        self.assertEqual(first.description, "Input VTK file")

    def test_positional_arg_invalid_kind_rejected(self) -> None:
        """A positional_arg with an invalid argument_kind must raise ValueError."""
        import tempfile
        extra = 'positional_args = [{name = "f", argument_kind = "integer", description = "x"}]\n'
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            with self.assertRaises(ValueError):
                _parse_manifest(p)


class TestProducesField(unittest.TestCase):
    """Phase 4: produces field on the manifest (Gap B coverage)."""

    def _make_produces_toml(self, produces_str: str) -> str:
        return produces_str

    def test_produces_well_formed_entry_accepted(self) -> None:
        """A well-formed produces entry must parse into a ProducesEntry tuple."""
        import tempfile
        # Use [[produces]] array-of-tables syntax (TOML 1.0 compatible).
        extra = textwrap.dedent(
            """\
            [[produces]]
            artifact_id = "graph_dict"
            path_pattern = "constant/purkinjeGraph"
            format = "openfoam_log"
            description = "Output dict"
            produced_by = "1DgraphToFoam"
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            manifest = _parse_manifest(p)
        self.assertEqual(len(manifest.produces), 1)
        entry = manifest.produces[0]
        self.assertIsInstance(entry, ProducesEntry)
        self.assertEqual(entry.artifact_id, "graph_dict")
        self.assertEqual(entry.path_pattern, "constant/purkinjeGraph")
        self.assertEqual(entry.format, "openfoam_log")
        self.assertFalse(entry.optional)
        self.assertFalse(entry.time_indexed)

    def test_produces_invalid_format_rejected(self) -> None:
        """A produces entry with an unknown format must raise ValueError."""
        import tempfile
        # Single-line inline table is TOML 1.0 compatible.
        extra = 'produces = [{artifact_id = "x", path_pattern = "out/file.txt", format = "xml_dump"}]\n'
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            with self.assertRaises(ValueError, msg="'xml_dump' should not be a valid format"):
                _parse_manifest(p)

    def test_produces_path_pattern_bad_placeholder_raises_at_load(self) -> None:
        """Gap B: a typo like {caseId} in path_pattern must raise ValueError at TOML-load time."""
        import tempfile
        # Use [[produces]] array-of-tables syntax (TOML 1.0 compatible).
        extra = textwrap.dedent(
            """\
            [[produces]]
            artifact_id = "x"
            path_pattern = "{caseId}/output.csv"
            format = "csv_probe"
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            with self.assertRaises(ValueError, msg="{caseId} should be caught at load time"):
                _parse_manifest(p)

    def test_produces_valid_placeholder_accepted(self) -> None:
        """path_pattern with {case_id} placeholder must parse without error."""
        import tempfile
        # Use [[produces]] array-of-tables syntax (TOML 1.0 compatible).
        extra = textwrap.dedent(
            """\
            [[produces]]
            artifact_id = "sweep_csv"
            path_pattern = "postProcessing/{case_id}_sweep.csv"
            format = "csv_sweep"
            time_indexed = false
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            manifest = _parse_manifest(p)
        self.assertEqual(manifest.produces[0].path_pattern, "postProcessing/{case_id}_sweep.csv")

    def test_produces_optional_fields_default(self) -> None:
        """optional and time_indexed default to False when not set."""
        import tempfile
        extra = 'produces = [{artifact_id="x", path_pattern="out.log", format="openfoam_log"}]\n'
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            manifest = _parse_manifest(p)
        e = manifest.produces[0]
        self.assertFalse(e.optional)
        self.assertFalse(e.time_indexed)
        self.assertEqual(e.variables, ())
        self.assertEqual(e.description, "")
        self.assertEqual(e.produced_by, "")


class TestOutputsProducesBackwardCompat(unittest.TestCase):
    """Phase 4: outputs stays valid; UserWarning fires on disagreement."""

    def test_outputs_alone_still_valid(self) -> None:
        """A manifest with only outputs (no produces) must load without warnings."""
        import tempfile
        extra = 'outputs = ["constant/purkinjeGraph"]\n'
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                manifest = _parse_manifest(p)
            compat_warnings = [x for x in w if issubclass(x.category, UserWarning)]
            self.assertEqual(len(compat_warnings), 0)
        self.assertEqual(manifest.outputs, ("constant/purkinjeGraph",))

    def test_outputs_and_produces_agreement_no_warning(self) -> None:
        """When outputs paths all appear in produces path_patterns, no warning is emitted."""
        import tempfile
        # Use [[produces]] array-of-tables syntax (TOML 1.0 compatible).
        extra = textwrap.dedent(
            """\
            outputs = ["constant/purkinjeGraph"]

            [[produces]]
            artifact_id = "graph_dict"
            path_pattern = "constant/purkinjeGraph"
            format = "openfoam_log"
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                _parse_manifest(p)
            compat_warnings = [x for x in w if issubclass(x.category, UserWarning)]
            self.assertEqual(len(compat_warnings), 0)

    def test_outputs_and_produces_disagreement_emits_warning(self) -> None:
        """When outputs has a path not in any produces entry, a UserWarning must be emitted."""
        import tempfile
        # Use [[produces]] array-of-tables syntax (TOML 1.0 compatible).
        extra = textwrap.dedent(
            """\
            outputs = ["constant/purkinjeGraph", "constant/otherDict"]

            [[produces]]
            artifact_id = "graph_dict"
            path_pattern = "constant/purkinjeGraph"
            format = "openfoam_log"
            """
        )
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "myutil"
            d.mkdir()
            p = d / "utility.manifest.toml"
            _write_minimal_toml(p, extra)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                _parse_manifest(p)
            compat_warnings = [x for x in w if issubclass(x.category, UserWarning)]
            self.assertGreater(
                len(compat_warnings),
                0,
                "Expected a UserWarning when outputs and produces disagree",
            )


class TestMigratedManifests(unittest.TestCase):
    """Phase 4 + 11a: all 9 manifests carry full schema (argument_kind + produces)."""

    # The three original Phase-4 manifests
    def test_1DgraphToFoam_produces_entries(self) -> None:
        """1DgraphToFoam manifest must have at least one produces entry after migration."""
        manifest = UTILITY_CATALOG.get("1DgraphToFoam")
        self.assertIsNotNone(manifest, "1DgraphToFoam not found in catalog")
        self.assertGreater(
            len(manifest.produces),
            0,
            "1DgraphToFoam must have at least one 'produces' entry after migration",
        )

    def test_1DgraphToFoam_flag_argument_kinds(self) -> None:
        """All 1DgraphToFoam flags must have a non-empty argument_kind after migration."""
        manifest = UTILITY_CATALOG.get("1DgraphToFoam")
        self.assertIsNotNone(manifest)
        for flag in manifest.flags:
            self.assertIn(
                flag.argument_kind,
                ALLOWED_ARGUMENT_KINDS,
                f"Flag {flag.name!r} has invalid argument_kind {flag.argument_kind!r}",
            )

    def test_runPurkinjeGraph_produces_entries(self) -> None:
        """runPurkinjeGraph manifest must have at least one produces entry."""
        manifest = UTILITY_CATALOG.get("runPurkinjeGraph")
        self.assertIsNotNone(manifest, "runPurkinjeGraph not found in catalog")
        self.assertGreater(len(manifest.produces), 0)

    def test_sweepCurrents_produces_entries(self) -> None:
        """sweepCurrents manifest must have at least one produces entry."""
        manifest = UTILITY_CATALOG.get("sweepCurrents")
        self.assertIsNotNone(manifest, "sweepCurrents not found in catalog")
        self.assertGreater(len(manifest.produces), 0)

    def test_migrated_produces_paths_match_outputs(self) -> None:
        """For the three original migrated utilities, every outputs path must appear in produces."""
        for util_name in ("1DgraphToFoam", "runPurkinjeGraph", "sweepCurrents"):
            manifest = UTILITY_CATALOG.get(util_name)
            self.assertIsNotNone(manifest, f"{util_name} not in catalog")
            produces_paths = {e.path_pattern for e in manifest.produces}
            for out_path in manifest.outputs:
                self.assertIn(
                    out_path,
                    produces_paths,
                    f"{util_name}: outputs path {out_path!r} not mirrored in produces",
                )

    # Phase 11a: six newly migrated manifests
    def test_checkMeshGeometry_produces_entries(self) -> None:
        """checkMeshGeometry must have at least one produces entry."""
        manifest = UTILITY_CATALOG.get("checkMeshGeometry")
        self.assertIsNotNone(manifest, "checkMeshGeometry not found in catalog")
        self.assertGreater(len(manifest.produces), 0)

    def test_checkMeshGeometry_flag_argument_kinds(self) -> None:
        """All checkMeshGeometry flags must carry a valid argument_kind."""
        manifest = UTILITY_CATALOG.get("checkMeshGeometry")
        self.assertIsNotNone(manifest)
        for flag in manifest.flags:
            self.assertIn(
                flag.argument_kind,
                ALLOWED_ARGUMENT_KINDS,
                f"Flag {flag.name!r} has invalid argument_kind {flag.argument_kind!r}",
            )

    def test_ionicHeterogeneityProbe_produces_entries(self) -> None:
        """ionicHeterogeneityProbe must have at least four produces entries (one per output CSV)."""
        manifest = UTILITY_CATALOG.get("ionicHeterogeneityProbe")
        self.assertIsNotNone(manifest, "ionicHeterogeneityProbe not found in catalog")
        self.assertGreaterEqual(len(manifest.produces), 4)

    def test_ionicHeterogeneityProbe_produces_paths_match_outputs(self) -> None:
        """Every ionicHeterogeneityProbe outputs path must appear in produces."""
        manifest = UTILITY_CATALOG.get("ionicHeterogeneityProbe")
        self.assertIsNotNone(manifest)
        produces_paths = {e.path_pattern for e in manifest.produces}
        for out_path in manifest.outputs:
            self.assertIn(
                out_path,
                produces_paths,
                f"ionicHeterogeneityProbe: outputs path {out_path!r} not mirrored in produces",
            )

    def test_listCellModelsVariables_produces_entries(self) -> None:
        """listCellModelsVariables must have at least one produces entry."""
        manifest = UTILITY_CATALOG.get("listCellModelsVariables")
        self.assertIsNotNone(manifest, "listCellModelsVariables not found in catalog")
        self.assertGreater(len(manifest.produces), 0)

    def test_listCellModelsVariables_produces_paths_match_outputs(self) -> None:
        """Every listCellModelsVariables outputs path must appear in produces."""
        manifest = UTILITY_CATALOG.get("listCellModelsVariables")
        self.assertIsNotNone(manifest)
        produces_paths = {e.path_pattern for e in manifest.produces}
        for out_path in manifest.outputs:
            self.assertIn(
                out_path,
                produces_paths,
                f"listCellModelsVariables: outputs path {out_path!r} not mirrored in produces",
            )

    def test_setTorsoOrganConductivityField_produces_entries(self) -> None:
        """setTorsoOrganConductivityField must have at least one produces entry."""
        manifest = UTILITY_CATALOG.get("setTorsoOrganConductivityField")
        self.assertIsNotNone(manifest, "setTorsoOrganConductivityField not found in catalog")
        self.assertGreater(len(manifest.produces), 0)

    def test_setTorsoOrganConductivityField_flag_argument_kinds(self) -> None:
        """All setTorsoOrganConductivityField flags must carry a valid argument_kind."""
        manifest = UTILITY_CATALOG.get("setTorsoOrganConductivityField")
        self.assertIsNotNone(manifest)
        for flag in manifest.flags:
            self.assertIn(
                flag.argument_kind,
                ALLOWED_ARGUMENT_KINDS,
                f"Flag {flag.name!r} has invalid argument_kind {flag.argument_kind!r}",
            )

    def test_setTorsoOrganConductivityField_produces_paths_match_outputs(self) -> None:
        """Every setTorsoOrganConductivityField outputs path must appear in produces."""
        manifest = UTILITY_CATALOG.get("setTorsoOrganConductivityField")
        self.assertIsNotNone(manifest)
        produces_paths = {e.path_pattern for e in manifest.produces}
        for out_path in manifest.outputs:
            self.assertIn(
                out_path,
                produces_paths,
                f"setTorsoOrganConductivityField: outputs path {out_path!r} not mirrored in produces",
            )

    def test_newVtkUnstructuredToFoam_produces_entries(self) -> None:
        """newVtkUnstructuredToFoam must have at least one produces entry."""
        manifest = UTILITY_CATALOG.get("newVtkUnstructuredToFoam")
        self.assertIsNotNone(manifest, "newVtkUnstructuredToFoam not found in catalog")
        self.assertGreater(len(manifest.produces), 0)

    def test_newVtkUnstructuredToFoam_flag_argument_kinds(self) -> None:
        """All newVtkUnstructuredToFoam flags must carry a valid argument_kind."""
        manifest = UTILITY_CATALOG.get("newVtkUnstructuredToFoam")
        self.assertIsNotNone(manifest)
        for flag in manifest.flags:
            self.assertIn(
                flag.argument_kind,
                ALLOWED_ARGUMENT_KINDS,
                f"Flag {flag.name!r} has invalid argument_kind {flag.argument_kind!r}",
            )

    def test_newVtkUnstructuredToFoam_positional_args(self) -> None:
        """newVtkUnstructuredToFoam must have at least one positional_arg (the vtk-file)."""
        manifest = UTILITY_CATALOG.get("newVtkUnstructuredToFoam")
        self.assertIsNotNone(manifest)
        self.assertGreater(len(manifest.positional_args), 0)
        self.assertEqual(manifest.positional_args[0].argument_kind, "path")

    def test_newVtkUnstructuredToFoam_produces_paths_match_outputs(self) -> None:
        """Every newVtkUnstructuredToFoam outputs path must appear in produces."""
        manifest = UTILITY_CATALOG.get("newVtkUnstructuredToFoam")
        self.assertIsNotNone(manifest)
        produces_paths = {e.path_pattern for e in manifest.produces}
        for out_path in manifest.outputs:
            self.assertIn(
                out_path,
                produces_paths,
                f"newVtkUnstructuredToFoam: outputs path {out_path!r} not mirrored in produces",
            )

    def test_setFibreField_produces_entries(self) -> None:
        """setFibreField must have produces entries covering all 11 outputs."""
        manifest = UTILITY_CATALOG.get("setFibreField")
        self.assertIsNotNone(manifest, "setFibreField not found in catalog")
        self.assertGreaterEqual(len(manifest.produces), 11)

    def test_setFibreField_produces_paths_match_outputs(self) -> None:
        """Every setFibreField outputs path must appear in produces."""
        manifest = UTILITY_CATALOG.get("setFibreField")
        self.assertIsNotNone(manifest)
        produces_paths = {e.path_pattern for e in manifest.produces}
        for out_path in manifest.outputs:
            self.assertIn(
                out_path,
                produces_paths,
                f"setFibreField: outputs path {out_path!r} not mirrored in produces",
            )

    def test_total_migrated_manifest_count(self) -> None:
        """All 9 utility manifests must carry at least one produces entry."""
        all_migrated = [
            "1DgraphToFoam",
            "runPurkinjeGraph",
            "sweepCurrents",
            "checkMeshGeometry",
            "ionicHeterogeneityProbe",
            "listCellModelsVariables",
            "setTorsoOrganConductivityField",
            "newVtkUnstructuredToFoam",
            "setFibreField",
        ]
        missing_produces = [
            name
            for name in all_migrated
            if UTILITY_CATALOG.get(name) is None
            or len(UTILITY_CATALOG[name].produces) == 0
        ]
        self.assertEqual(
            missing_produces,
            [],
            f"These manifests are not yet fully migrated (missing produces): {missing_produces}",
        )


if __name__ == "__main__":
    unittest.main()
