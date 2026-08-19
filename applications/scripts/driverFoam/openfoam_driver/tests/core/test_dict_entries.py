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
#     test_dict_entries
#
# Description
#     Tests dict entries logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.dict_entries import (
    get_electro_property_entry_groups,
    PHYSICS_PROPERTY_ENTRIES,
    all_documented_driver_paths,
)
from openfoam_driver.plugins.cardiacfoam.overrides import apply_electro_property_overrides
from openfoam_driver.tests.conftest import assert_foam_entry


class TestDictEntryCatalog(unittest.TestCase):
    def test_catalog_contains_core_physics_and_electro_paths(self) -> None:
        physics_paths = {entry.driver_path for entry in PHYSICS_PROPERTY_ENTRIES}
        self.assertEqual(physics_paths, {"type"})

        documented = set(all_documented_driver_paths())
        expected = {
            "myocardiumSolver",
            "$ELECTRO_MODEL_COEFFS.solutionAlgorithm",
            "$ELECTRO_MODEL_COEFFS.ionicModel",
            "$ELECTRO_MODEL_COEFFS.tissue",
            "$ELECTRO_MODEL_COEFFS.dimension",
            "$ELECTRO_MODEL_COEFFS.writeAfterTime",
            "$ELECTRO_MODEL_COEFFS.utilities",
            "$ELECTRO_MODEL_COEFFS.initSampleCell",
            "$ELECTRO_MODEL_COEFFS.outputVariables.ionic.export",
            "$ELECTRO_MODEL_COEFFS.outputVariables.activeTension.export",
            "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S1",
            "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensity",
            "$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach",
            "$ELECTRO_MODEL_COEFFS.bathPotentialDomain.bathCellZones",
            "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver",
            "$ELECTRO_MODEL_COEFFS.activeTensionModel",
            "$ELECTRO_MODEL_COEFFS.couplingSignal",
        }
        self.assertTrue(expected.issubset(documented))

    def test_catalog_paths_are_unique(self) -> None:
        documented = all_documented_driver_paths()
        self.assertEqual(len(documented), len(set(documented)))

    def test_catalog_mentions_existing_source_files(self) -> None:
        import pytest
        repo_root = Path(__file__).resolve()
        _found = None
        for parent in repo_root.parents:
            if (parent / "src").exists() and (parent / "applications").exists():
                _found = parent
                break
        if _found is None:
            pytest.skip(
                "Monorepo src/ tree not present — skipping source-ref file-existence check. "
                "Run from the full cardiacFoam checkout to enable this test."
            )
        repo_root = _found

        for entry in PHYSICS_PROPERTY_ENTRIES:
            for source_ref in entry.source_refs:
                self.assertTrue((repo_root / source_ref).exists(), source_ref)

        for entries in get_electro_property_entry_groups().values():
            for entry in entries:
                for source_ref in entry.source_refs:
                    self.assertTrue((repo_root / source_ref).exists(), source_ref)

    def test_catalog_exposes_gui_value_hints_for_key_entries(self) -> None:
        type_entry = PHYSICS_PROPERTY_ENTRIES[0]
        self.assertEqual(type_entry.value_kind, "enum")
        self.assertIn("electroMechanicalModel", type_entry.enum_values)

        monodomain_entries = {
            entry.driver_path: entry for entry in get_electro_property_entry_groups()["monodomain"]
        }
        self.assertEqual(
            monodomain_entries["$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMin"].value_kind,
            "vector3",
        )
        self.assertEqual(
            monodomain_entries["$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensity"].value_kind,
            "dimensioned_scalar_literal",
        )

        ecg_entries = {entry.driver_path: entry for entry in get_electro_property_entry_groups()["ecg"]}
        self.assertTrue(
            ecg_entries[
                "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.electrodePositions.<electrode>"
            ].dynamic_path
        )


class TestDeepElectroOverrides(unittest.TestCase):
    def test_apply_electro_property_overrides_updates_dimensioned_and_dynamic_entries(self) -> None:
        text = "\n".join(
            [
                "myocardiumSolver monodomainSolver;",
                "",
                "monodomainSolverCoeffs",
                "{",
                "    conductivity [-1 -3 3 0 0 2 0] (0.133 0 0 0.017 0 0.017);",
                "    externalStimulus",
                "    {",
                "        stimulusIntensity [0 -3 0 0 0 1 0] 50000;",
                "    }",
                "    ecgDomains",
                "    {",
                "        ECG",
                "        {",
                "            ecgSolver pseudoECG;",
                "            electrodePositions",
                "                {",
                    "                    V1 (-0.02 -0.28 -0.07);",
                "                }",
                "        }",
                "    }",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            apply_electro_property_overrides(
                path,
                {
                    "$ELECTRO_MODEL_COEFFS.conductivity": "[-1 -3 3 0 0 2 0] (0.2 0 0 0.03 0 0.03)",
                    "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensity": "[0 -3 0 0 0 1 0] 75000",
                    "$ELECTRO_MODEL_COEFFS.ecgDomains.ECG.electrodePositions.V1": "(1 2 3)",
                },
            )

            coeffs = "monodomainSolverCoeffs"
            assert_foam_entry(
                path,
                "conductivity",
                "[-1 -3 3 0 0 2 0] (0.2 0 0 0.03 0 0.03)",
                scope=coeffs,
            )
            assert_foam_entry(
                path,
                "stimulusIntensity",
                "[0 -3 0 0 0 1 0] 75000",
                scope=(coeffs, "externalStimulus"),
            )
            assert_foam_entry(
                path,
                "V1",
                "(1 2 3)",
                scope=(coeffs, "ecgDomains", "ECG", "electrodePositions"),
            )


class TestConductionSystemSchemaContract(unittest.TestCase):
    """Verifies that the conduction_system group uses the keys the C++ code actually reads."""

    def setUp(self):
        self.entries = {
            e.driver_path: e
            for e in get_electro_property_entry_groups()["conduction_system"]
        }

    def test_conduction_domain_selector_key_is_conductionSystemDomain(self):
        # C++ uses lowercase key names in dictionary lookups.
        matching = [
            p for p in self.entries
            if p.endswith(".conductionSystemDomain")
        ]
        self.assertTrue(
            len(matching) >= 1,
            "Expected at least one entry whose path ends with '.conductionSystemDomain'"
        )

    def test_graph_file_schema_is_documented_on_the_coeffs_subdict(self):
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.graphFile",
            self.entries,
        )
        self.assertNotIn(
            "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.pvjNodes",
            self.entries,
        )
        self.assertNotIn(
            "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.pvjLocations",
            self.entries,
        )

    def test_root_stimulus_sub_entries_documented(self):
        for sub in ("startTime", "duration", "intensity", "node"):
            matching = [
                p
                for p in self.entries
                if p.endswith(f".purkinjeGraphModelCoeffs.rootStimulus.{sub}")
            ]
            self.assertTrue(
                len(matching) >= 1,
                f"rootStimulus.{sub} not documented"
            )

    def test_purkinjeGraphModelCoeffs_chi_and_cm_documented(self):
        chi_keys = [p for p in self.entries if p.endswith(".purkinjeGraphModelCoeffs.chi")]
        cm_keys  = [p for p in self.entries if p.endswith(".purkinjeGraphModelCoeffs.cm")]
        self.assertTrue(len(chi_keys) >= 1, "purkinjeGraphModelCoeffs.chi not documented")
        self.assertTrue(len(cm_keys)  >= 1, "purkinjeGraphModelCoeffs.cm not documented")

class TestDomainCouplingSchemaContract(unittest.TestCase):
    """Verifies that the domain_couplings group owns the domainCouplings schema."""

    def setUp(self):
        self.entries = {
            e.driver_path: e
            for e in get_electro_property_entry_groups()["domain_couplings"]
        }

    def test_coupler_selector_key_is_electroDomainCoupler(self):
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.electroDomainCoupler",
            self.entries,
        )

    def test_coupling_helper_keys_documented(self):
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.conductionNetworkDomain",
            self.entries,
        )
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.rPvj",
            self.entries,
        )
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjRadius",
            self.entries,
        )
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjCouplingScheme",
            self.entries,
        )
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.couplingMode",
            self.entries,
        )

    def test_common_model_coeffs_owns_electrophysics_advance_scheme(self):
        common_entries = {
            e.driver_path: e
            for e in get_electro_property_entry_groups()["common_model_coeffs"]
        }
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.electrophysicsAdvanceScheme",
            common_entries,
        )
        self.assertEqual(
            common_entries["$ELECTRO_MODEL_COEFFS.electrophysicsAdvanceScheme"].enum_values,
            ("staggeredElectrophysicsAdvanceScheme",),
        )


class TestDictEntryStructuredConstraints(unittest.TestCase):
    """DictEntry exposes four structured-constraint fields so that
    constraints can be expressed in a form the validator can evaluate.

    The fields are additive (P8 additive-only policy): every existing
    DictEntry must construct unchanged with empty defaults.
    """

    def _build_entry(self, **overrides) -> "DictEntry":  # noqa: F821 - imported below
        from openfoam_driver.dict_entries import DictEntry
        defaults = {
            "driver_path": "test.path",
            "description": "fixture",
            "source_refs": ("ref.C",),
        }
        defaults.update(overrides)
        return DictEntry(**defaults)

    def test_applicable_when_defaults_empty(self) -> None:
        entry = self._build_entry()
        self.assertEqual(entry.applicable_when, {})

    def test_forbidden_when_defaults_empty(self) -> None:
        entry = self._build_entry()
        self.assertEqual(entry.forbidden_when, {})

    def test_required_when_defaults_empty(self) -> None:
        entry = self._build_entry()
        self.assertEqual(entry.required_when, {})

    def test_mutually_exclusive_with_defaults_empty(self) -> None:
        entry = self._build_entry()
        self.assertEqual(entry.mutually_exclusive_with, ())

    def test_applicable_when_accepts_value_predicate(self) -> None:
        entry = self._build_entry(
            applicable_when={"myocardiumSolver": "monodomainSolver"},
        )
        self.assertEqual(
            entry.applicable_when, {"myocardiumSolver": "monodomainSolver"},
        )

    def test_applicable_when_accepts_value_list_predicate(self) -> None:
        """Some constraints target multiple legal values
        (e.g. 'manufactured ionic models X, Y, Z')."""
        entry = self._build_entry(
            applicable_when={
                "ionicModel": (
                    "monodomainFDAManufactured",
                    "bidomainFDAManufactured",
                    "bathBidomainFDAManufactured",
                ),
            },
        )
        self.assertEqual(len(entry.applicable_when["ionicModel"]), 3)

    def test_forbidden_when_accepts_value_predicate(self) -> None:
        entry = self._build_entry(
            forbidden_when={"myocardiumSolver": "eikonalSolver"},
        )
        self.assertEqual(entry.forbidden_when["myocardiumSolver"], "eikonalSolver")

    def test_required_when_accepts_value_predicate(self) -> None:
        entry = self._build_entry(
            required_when={"myocardiumSolver": "singleCellSolver"},
        )
        self.assertEqual(
            entry.required_when["myocardiumSolver"], "singleCellSolver",
        )

    def test_mutually_exclusive_with_accepts_path_tuple(self) -> None:
        entry = self._build_entry(
            mutually_exclusive_with=("stimulusDurationList",),
        )
        self.assertEqual(entry.mutually_exclusive_with, ("stimulusDurationList",))

    def test_entry_remains_frozen(self) -> None:
        """The additive fields must not loosen the existing
        immutability guarantee on DictEntry."""
        import dataclasses
        entry = self._build_entry()
        with self.assertRaises(dataclasses.FrozenInstanceError):
            entry.applicable_when = {"x": "y"}  # type: ignore[misc]

    def test_existing_entries_in_catalog_have_empty_defaults(self) -> None:
        """Every entry in the live catalog must still construct cleanly
        with empty structured-constraint fields — migration is opt-in
        per-entry, not a forced rewrite."""
        from openfoam_driver.dict_entries import (
            get_electro_property_entry_groups,
            PHYSICS_PROPERTY_ENTRIES,
        )
        all_entries = list(PHYSICS_PROPERTY_ENTRIES)
        for group in get_electro_property_entry_groups().values():
            all_entries.extend(group)
        self.assertGreater(len(all_entries), 80)  # sanity: we have 87+ today
        for entry in all_entries:
            # No AttributeError accessing the new fields.
            self.assertIsInstance(entry.applicable_when, dict)
            self.assertIsInstance(entry.forbidden_when, dict)
            self.assertIsInstance(entry.required_when, dict)
            self.assertIsInstance(entry.mutually_exclusive_with, tuple)


class TestElectroPropertiesPresenceScans(unittest.TestCase):
    """Tests for the three presence helpers used by the predictor's
    domain-aware handlers.
    """

    def _write(self, body: str) -> Path:
        import tempfile
        from pathlib import Path
        temp = tempfile.mkdtemp()
        path = Path(temp) / "electroProperties"
        path.write_text(body)
        return path

    def test_has_block_finds_top_level_block(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import electro_properties_has_block
        path = self._write(
            "myocardiumSolver bidomainSolver;\n"
            "bidomainSolverCoeffs\n{\n  ionicModel TNNP;\n}\n"
            "ecgDomains\n{\n  myECG { ecgSolver pseudoECG; }\n}\n"
        )
        self.assertTrue(electro_properties_has_block(path, "ecgDomains"))
        self.assertFalse(electro_properties_has_block(path, "conductionNetworkDomains"))

    def test_has_block_handles_inline_brace(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import electro_properties_has_block
        path = self._write(
            "myocardiumSolver monodomainSolver;\n"
            "conductionNetworkDomains { purk { } }\n"
        )
        self.assertTrue(
            electro_properties_has_block(path, "conductionNetworkDomains")
        )

    def test_has_block_ignores_substring_matches(self) -> None:
        """The scan must match block declarations, not keys whose names
        happen to contain the target word."""
        from openfoam_driver.plugins.cardiacfoam.detection import electro_properties_has_block
        path = self._write(
            "myocardiumSolver monodomainSolver;\n"
            "monodomainSolverCoeffs\n{\n  ecgDomainsCount 0;\n}\n"
        )
        self.assertFalse(electro_properties_has_block(path, "ecgDomains"))

    def test_detect_verification_model_type_present(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import detect_verification_model_type
        path = self._write(
            "myocardiumSolver monodomainSolver;\n"
            "monodomainSolverCoeffs\n{\n"
            "  ionicModel monodomainFDAManufactured;\n"
            "  verificationModel\n  {\n"
            "    type manufacturedFDAMonodomainVerifier;\n"
            "  }\n"
            "}\n"
        )
        self.assertEqual(
            detect_verification_model_type(path),
            "manufacturedFDAMonodomainVerifier",
        )

    def test_detect_verification_model_type_absent_returns_none(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import detect_verification_model_type
        path = self._write(
            "myocardiumSolver monodomainSolver;\n"
            "monodomainSolverCoeffs\n{\n  ionicModel TNNP;\n}\n"
        )
        self.assertIsNone(detect_verification_model_type(path))


class TestDetectActiveTensionModelName(unittest.TestCase):
    def _write(self, text: str) -> Path:
        p = Path(tempfile.mkdtemp()) / "electroProperties"
        p.write_text(text)
        return p

    def test_detects_nash_panfilov(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import detect_active_tension_model_name
        props = self._write(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n{\n"
            "    activeTensionModel NashPanfilov;\n"
            "}\n"
        )
        self.assertEqual(detect_active_tension_model_name(props), "NashPanfilov")

    def test_detects_goktepe_kuhl(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import detect_active_tension_model_name
        props = self._write(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n{\n"
            "    activeTensionModel GoktepeKuhl;\n"
            "}\n"
        )
        self.assertEqual(detect_active_tension_model_name(props), "GoktepeKuhl")

    def test_returns_none_when_block_absent(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import detect_active_tension_model_name
        props = self._write(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n{\n"
            "    ionicModel TNNP;\n"
            "}\n"
        )
        self.assertIsNone(detect_active_tension_model_name(props))


class TestDetectActiveTensionExportList(unittest.TestCase):
    def _write(self, text: str) -> Path:
        p = Path(tempfile.mkdtemp()) / "electroProperties"
        p.write_text(text)
        return p

    def test_detects_ta_export(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import detect_active_tension_export_list
        props = self._write(
            "myocardiumSolver monodomainSolver;\n"
            "monodomainSolverCoeffs\n{\n"
            "    outputVariables\n    {\n"
            "        activeTension\n        {\n"
            "            export ( Ta );\n"
            "        }\n"
            "    }\n"
            "}\n"
        )
        self.assertEqual(detect_active_tension_export_list(props), ("Ta",))

    def test_returns_none_when_absent(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import detect_active_tension_export_list
        props = self._write(
            "myocardiumSolver monodomainSolver;\n"
            "monodomainSolverCoeffs\n{\n"
            "    ionicModel TNNP;\n"
            "}\n"
        )
        self.assertIsNone(detect_active_tension_export_list(props))


class TestEmptyExportListIsKnownEmpty(unittest.TestCase):
    """An explicit ``export ()`` means "known, and empty" -- never "unknown".

    Regression guard for the artifact-prediction bug in which an empty token
    tuple was coerced to ``None`` (``tokens if tokens else None``).  ``None``
    is the predictor's signal that nothing was declared, so it fell back to
    the ionic/active-tension catalog's recommended exports and predicted
    artifacts the solver was explicitly told not to write, failing the strict
    run with ``missing_expected_artifacts``.  The detector's ``None`` must
    mean "no ``export`` block found" and nothing else.
    """

    def _case_root(self, text: str) -> Path:
        case_root = Path(tempfile.mkdtemp())
        constant = case_root / "constant"
        constant.mkdir()
        (constant / "electroProperties").write_text(text)
        return case_root

    def _electro_properties(self, *, ionic_export: str, at_export: str | None = None) -> str:
        at_block = ""
        if at_export is not None:
            at_block = (
                "        activeTension\n        {\n"
                f"            export ( {at_export} );\n"
                "        }\n"
            )
        return (
            "myocardiumSolver monodomainSolver;\n"
            "monodomainSolverCoeffs\n{\n"
            "    ionicModel AlievPanfilov;\n"
            "    outputVariables\n    {\n"
            "        ionic\n        {\n"
            f"            export ( {ionic_export} );\n"
            "        }\n"
            f"{at_block}"
            "    }\n"
            "}\n"
        )

    def test_empty_ionic_export_detects_as_empty_tuple_not_none(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import detect_ionic_export_list

        case_root = self._case_root(self._electro_properties(ionic_export=""))
        declared = detect_ionic_export_list(case_root / "constant" / "electroProperties")
        self.assertIsNotNone(declared, "empty export () must not be reported as unknown")
        self.assertEqual(declared, ())

    def test_empty_active_tension_export_detects_as_empty_tuple_not_none(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.detection import (
            detect_active_tension_export_list,
        )

        case_root = self._case_root(
            self._electro_properties(ionic_export="Vm", at_export="")
        )
        declared = detect_active_tension_export_list(
            case_root / "constant" / "electroProperties"
        )
        self.assertIsNotNone(declared, "empty export () must not be reported as unknown")
        self.assertEqual(declared, ())

    def test_predictor_honours_empty_export_over_catalog_defaults(self) -> None:
        """The behaviour the detector fix exists to protect."""
        from openfoam_driver.plugins.cardiacfoam.artifacts_predictor import (
            _exported_ionic_variables,
        )
        from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import (
            IONIC_MODEL_CATALOG,
        )

        # The catalog must actually recommend something, or this test would
        # pass for the wrong reason.
        self.assertTrue(IONIC_MODEL_CATALOG["AlievPanfilov"].recommended_exports)

        empty = self._case_root(self._electro_properties(ionic_export=""))
        self.assertEqual(_exported_ionic_variables(empty, "AlievPanfilov"), ())

    def test_predictor_still_falls_back_when_export_block_absent(self) -> None:
        """The other half of the distinction: absent really is unknown."""
        from openfoam_driver.plugins.cardiacfoam.artifacts_predictor import (
            _exported_ionic_variables,
        )
        from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import (
            IONIC_MODEL_CATALOG,
        )

        case_root = self._case_root(
            "myocardiumSolver monodomainSolver;\n"
            "monodomainSolverCoeffs\n{\n"
            "    ionicModel AlievPanfilov;\n"
            "}\n"
        )
        self.assertEqual(
            _exported_ionic_variables(case_root, "AlievPanfilov"),
            IONIC_MODEL_CATALOG["AlievPanfilov"].recommended_exports,
        )


class TestReadFoamEntry(unittest.TestCase):
    """read_foam_entry returns the raw value string for a key, or None."""

    def test_reads_top_level_key(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f"
            p.write_text("myocardiumSolver monodomainSolver;\n")
            from openfoam_driver.core.runtime.mutators import read_foam_entry
            result = read_foam_entry(p, "myocardiumSolver")
            self.assertEqual(result, "monodomainSolver")

    def test_reads_scoped_key(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f"
            p.write_text(
                "myocardiumSolver monodomainSolver;\n"
                "monodomainSolverCoeffs\n"
                "{\n"
                "    ionicModel TNNP;\n"
                "    solutionAlgorithm implicit;\n"
                "}\n"
            )
            from openfoam_driver.core.runtime.mutators import read_foam_entry
            result = read_foam_entry(p, "ionicModel", scope="monodomainSolverCoeffs")
            self.assertEqual(result, "TNNP")

    def test_reads_nested_scope_key(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f"
            p.write_text(
                "singleCellSolverCoeffs\n"
                "{\n"
                "    singleCellStimulus\n"
                "    {\n"
                "        stim_amplitude 60;\n"
                "    }\n"
                "}\n"
            )
            from openfoam_driver.core.runtime.mutators import read_foam_entry
            result = read_foam_entry(
                p, "stim_amplitude",
                scope=["singleCellSolverCoeffs", "singleCellStimulus"],
            )
            self.assertEqual(result, "60")

    def test_missing_key_returns_none(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f"
            p.write_text("myocardiumSolver monodomainSolver;\n")
            from openfoam_driver.core.runtime.mutators import read_foam_entry
            result = read_foam_entry(p, "notAKey")
            self.assertIsNone(result)

    def test_missing_scope_returns_none(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f"
            p.write_text("myocardiumSolver monodomainSolver;\n")
            from openfoam_driver.core.runtime.mutators import read_foam_entry
            result = read_foam_entry(p, "ionicModel", scope="noSuchBlock")
            self.assertIsNone(result)

    def test_strips_inline_comment_from_value(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f"
            p.write_text("solutionAlgorithm implicit; // time discretisation\n")
            from openfoam_driver.core.runtime.mutators import read_foam_entry
            result = read_foam_entry(p, "solutionAlgorithm")
            self.assertEqual(result, "implicit")

    def test_nonexistent_file_returns_none(self) -> None:
        from openfoam_driver.core.runtime.mutators import read_foam_entry
        result = read_foam_entry(Path("/no/such/file"), "key")
        self.assertIsNone(result)


class TestUpdateControlDict(unittest.TestCase):
    """update_control_dict patches deltaT / endTime in an existing controlDict."""

    @staticmethod
    def _write_control_dict(d: str) -> Path:
        p = Path(d) / "controlDict"
        p.write_text("deltaT    0.05;\nendTime   1.0;\n")
        return p

    def test_updates_delta_t(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = self._write_control_dict(d)
            from openfoam_driver.core.runtime.mutators import update_control_dict
            update_control_dict(p, delta_t=0.001)
            self.assertIn("0.001", p.read_text())

    def test_updates_end_time(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = self._write_control_dict(d)
            from openfoam_driver.core.runtime.mutators import update_control_dict
            update_control_dict(p, end_time=0.5)
            self.assertIn("0.5", p.read_text())

    def test_updates_both_in_one_call(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = self._write_control_dict(d)
            from openfoam_driver.core.runtime.mutators import update_control_dict
            update_control_dict(p, delta_t=0.002, end_time=0.1)
            text = p.read_text()
            self.assertIn("0.002", text)
            self.assertIn("0.1", text)

    def test_none_args_leave_file_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = self._write_control_dict(d)
            original = p.read_text()
            from openfoam_driver.core.runtime.mutators import update_control_dict
            update_control_dict(p)
            self.assertEqual(p.read_text(), original)

    def test_missing_file_raises_file_not_found(self) -> None:
        from openfoam_driver.core.runtime.mutators import update_control_dict
        with self.assertRaises(FileNotFoundError):
            update_control_dict(Path("/no/such/controlDict"), delta_t=0.001)


class TestControlDictEntries(unittest.TestCase):
    """CONTROL_DICT_ENTRIES catalog shape contract."""

    def test_catalog_exposes_delta_t_and_end_time(self) -> None:
        from openfoam_driver.dict_entries import CONTROL_DICT_ENTRIES
        driver_paths = {e.driver_path for e in CONTROL_DICT_ENTRIES}
        self.assertIn("deltaT", driver_paths)
        self.assertIn("endTime", driver_paths)

    def test_time_entries_carry_seconds_unit(self) -> None:
        from openfoam_driver.dict_entries import CONTROL_DICT_ENTRIES
        time_entries = {"deltaT", "endTime", "startTime", "writeInterval"}
        for entry in CONTROL_DICT_ENTRIES:
            if entry.driver_path in time_entries:
                self.assertTrue(
                    entry.unit.startswith("s"),
                    f"{entry.driver_path} must carry a seconds unit, got '{entry.unit}'"
                )

    def test_entries_belong_to_solver_phase(self) -> None:
        from openfoam_driver.dict_entries import CONTROL_DICT_ENTRIES
        for entry in CONTROL_DICT_ENTRIES:
            self.assertIn("solver", entry.phases,
                          f"{entry.driver_path} must be in solver phase")

    def test_entries_are_marked_required(self) -> None:
        from openfoam_driver.dict_entries import CONTROL_DICT_ENTRIES
        for entry in CONTROL_DICT_ENTRIES:
            self.assertTrue(entry.required,
                            f"{entry.driver_path} must be required=True")


if __name__ == "__main__":
    unittest.main()


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
#     test_dict_entries_coverage
#
# Description
#     Tests dict entries coverage logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Coverage tests for ``DictEntry.phases`` classification.

Enforces that every catalogued entry declares at least one workflow phase, and
that every declared value is a valid ``Phase`` literal.
"""


import typing

from openfoam_driver.dict_entries import (
    get_electro_property_entry_groups,
    PHYSICS_PROPERTY_ENTRIES,
    Phase,
)

VALID_PHASES = set(typing.get_args(Phase))


def _all_entries():
    yield from PHYSICS_PROPERTY_ENTRIES
    for group in get_electro_property_entry_groups().values():
        yield from group


def test_every_dict_entry_has_at_least_one_phase():
    unclassified = [e for e in _all_entries() if not e.phases]
    assert not unclassified, (
        f"{len(unclassified)} entries have no phases: "
        + ", ".join(e.driver_path for e in unclassified[:10])
    )


def test_every_phase_value_is_a_valid_literal():
    invalid = []
    for e in _all_entries():
        bad = [p for p in e.phases if p not in VALID_PHASES]
        if bad:
            invalid.append((e.driver_path, bad))
    assert not invalid, f"invalid phases: {invalid}"

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
#     test_dict_entry_phases
#
# Description
#     Tests dict entry phases logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the ``Phase`` literal and multi-phase ``DictEntry.phases`` field.

Covers the DictEntry dataclass contract that lets a single entry belong to more
than one workflow phase.
"""


import typing

from openfoam_driver.dict_entries import DictEntry, Phase

VALID_PHASES = {"anatomy", "physics", "stimulus", "solver"}


def test_dict_entry_has_phases_field_accepting_a_frozenset():
    entry = DictEntry(
        driver_path="foo",
        description="x",
        source_refs=("bar",),
        phases=frozenset({"physics"}),
    )
    assert entry.phases == frozenset({"physics"})


def test_dict_entry_phases_supports_multi_phase_ownership():
    entry = DictEntry(
        driver_path="nRegions",
        description="number of regions",
        source_refs=("bar",),
        phases=frozenset({"anatomy", "solver"}),
    )
    assert entry.phases == frozenset({"anatomy", "solver"})


def test_dict_entry_phases_default_is_empty_frozenset():
    entry = DictEntry(driver_path="foo", description="x", source_refs=("bar",))
    assert entry.phases == frozenset()


def test_phase_literal_values():
    args = typing.get_args(Phase)
    assert set(args) == VALID_PHASES
