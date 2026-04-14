from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.dict_entries import (
    ELECTRO_PROPERTY_ENTRY_GROUPS,
    PHYSICS_PROPERTY_ENTRIES,
    all_documented_driver_paths,
)
from openfoam_driver.specs.common import apply_electro_property_overrides


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
            "$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusIntensity",
            "$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach",
            "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver",
            "$ELECTRO_MODEL_COEFFS.activeTensionModel.activeTensionModel",
            "$ELECTRO_MODEL_COEFFS.activeTensionModel.couplingSignal",
        }
        self.assertTrue(expected.issubset(documented))

    def test_catalog_paths_are_unique(self) -> None:
        documented = all_documented_driver_paths()
        self.assertEqual(len(documented), len(set(documented)))

    def test_catalog_mentions_existing_source_files(self) -> None:
        repo_root = Path(__file__).resolve()
        for parent in repo_root.parents:
            if (parent / "src").exists() and (parent / "applications").exists():
                repo_root = parent
                break
        else:
            self.fail("Could not locate repository root from test path")

        for entry in PHYSICS_PROPERTY_ENTRIES:
            for source_ref in entry.source_refs:
                self.assertTrue((repo_root / source_ref).exists(), source_ref)

        for entries in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
            for entry in entries:
                for source_ref in entry.source_refs:
                    self.assertTrue((repo_root / source_ref).exists(), source_ref)

    def test_catalog_exposes_gui_value_hints_for_key_entries(self) -> None:
        type_entry = PHYSICS_PROPERTY_ENTRIES[0]
        self.assertEqual(type_entry.value_kind, "enum")
        self.assertEqual(type_entry.ui_control, "select")
        self.assertIn("electroMechanicalModel", type_entry.enum_values)

        monodomain_entries = {
            entry.driver_path: entry for entry in ELECTRO_PROPERTY_ENTRY_GROUPS["monodomain"]
        }
        self.assertEqual(
            monodomain_entries["$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusLocationMin"].value_kind,
            "vector3",
        )
        self.assertEqual(
            monodomain_entries["$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusIntensity"].value_kind,
            "dimensioned_scalar_literal",
        )

        ecg_entries = {entry.driver_path: entry for entry in ELECTRO_PROPERTY_ENTRY_GROUPS["ecg"]}
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
                "    monodomainStimulus",
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
                    "$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusIntensity": "[0 -3 0 0 0 1 0] 75000",
                    "$ELECTRO_MODEL_COEFFS.ecgDomains.ECG.electrodePositions.V1": "(1 2 3)",
                },
            )

            updated = path.read_text()
            self.assertIn(
                "conductivity    [-1 -3 3 0 0 2 0] (0.2 0 0 0.03 0 0.03);",
                updated,
            )
            self.assertIn(
                "stimulusIntensity    [0 -3 0 0 0 1 0] 75000;",
                updated,
            )
            self.assertIn("V1    (1 2 3);", updated)


class TestConductionSystemSchemaContract(unittest.TestCase):
    """Verifies that the conduction_system group uses the keys the C++ code actually reads."""

    def setUp(self):
        self.entries = {
            e.driver_path: e
            for e in ELECTRO_PROPERTY_ENTRY_GROUPS["conduction_system"]
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

    def test_coupler_selector_key_is_electroDomainCoupler(self):
        # C++ uses lowercase key names in dictionary lookups.
        matching = [
            p for p in self.entries
            if p.endswith(".electroDomainCoupler")
        ]
        self.assertTrue(
            len(matching) >= 1,
            "Expected at least one entry whose path ends with '.electroDomainCoupler'"
        )

    def test_advance_scheme_key_is_electrophysicsAdvanceScheme(self):
        # C++ electrophysicsAdvanceScheme.C: dict.lookupOrDefault<word>("electrophysicsAdvanceScheme", ...)
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.electrophysicsAdvanceScheme",
            self.entries,
        )

    def test_advance_scheme_enum_includes_pimple_staggered(self):
        entry = self.entries["$ELECTRO_MODEL_COEFFS.electrophysicsAdvanceScheme"]
        self.assertIn(
            "pimpleStaggeredElectrophysicsAdvanceScheme",
            entry.enum_values,
        )

    def test_pvjNodes_and_pvjLocations_documented(self):
        pvj_nodes_keys = [p for p in self.entries if p.endswith(".pvjNodes")]
        pvj_locs_keys  = [p for p in self.entries if p.endswith(".pvjLocations")]
        self.assertTrue(len(pvj_nodes_keys) >= 1, "pvjNodes not documented")
        self.assertTrue(len(pvj_locs_keys)  >= 1, "pvjLocations not documented")

    def test_root_stimulus_sub_entries_documented(self):
        for sub in ("startTime", "duration", "intensity"):
            matching = [p for p in self.entries if p.endswith(f".rootStimulus.{sub}")]
            self.assertTrue(
                len(matching) >= 1,
                f"rootStimulus.{sub} not documented"
            )

    def test_purkinjeNetworkModelCoeffs_chi_and_cm_documented(self):
        chi_keys = [p for p in self.entries if p.endswith(".purkinjeNetworkModelCoeffs.chi")]
        cm_keys  = [p for p in self.entries if p.endswith(".purkinjeNetworkModelCoeffs.cm")]
        self.assertTrue(len(chi_keys) >= 1, "purkinjeNetworkModelCoeffs.chi not documented")
        self.assertTrue(len(cm_keys)  >= 1, "purkinjeNetworkModelCoeffs.cm not documented")

    def test_coupling_helper_keys_documented(self):
        network_domain_keys = [p for p in self.entries if p.endswith(".conductionNetworkDomain")]
        self.assertTrue(len(network_domain_keys) >= 1, "conductionNetworkDomain not documented")


if __name__ == "__main__":
    unittest.main()
