"""Tests for the dict_builder (plan §9.1).

The builder synthesizes a complete electroProperties dict from
minimum-viable agent intent (selectors + overrides). It enforces the same
constraints the validator does, so its output is guaranteed validator-clean
by construction.
"""
from __future__ import annotations

import unittest


class TestDictBuilderModule(unittest.TestCase):
    """Module-level structural contract — the import path + signature."""

    def test_module_exposes_build_electro_properties(self) -> None:
        from openfoam_driver.specs.dict_builder import build_electro_properties
        self.assertTrue(callable(build_electro_properties))

    def test_function_accepts_documented_kwargs(self) -> None:
        import inspect
        from openfoam_driver.specs.dict_builder import build_electro_properties
        sig = inspect.signature(build_electro_properties)
        params = sig.parameters
        self.assertIn("selectors", params)
        self.assertIn("overrides", params)
        self.assertIn("typical_value_fallback", params)
        # selectors is required positional/keyword; overrides + fallback are keyword-only.
        self.assertEqual(params["overrides"].kind, inspect.Parameter.KEYWORD_ONLY)
        self.assertEqual(params["typical_value_fallback"].kind, inspect.Parameter.KEYWORD_ONLY)


class TestMinimalSingleCellBuild(unittest.TestCase):
    """First behavioural test: a singleCellSolver + AlievPanfilov case must
    produce a string containing the FoamFile preamble and the chosen solver
    selector. Drives the bare-minimum end-to-end pipeline."""

    def test_returns_string_with_foamfile_preamble(self) -> None:
        from openfoam_driver.specs.dict_builder import build_electro_properties
        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "singleCellSolver",
                "ionicModel": "AlievPanfilov",
                "tissue": "myocyte",
            },
        )
        self.assertIsInstance(text, str)
        self.assertIn("FoamFile", text)
        self.assertIn("myocardiumSolver", text)
        self.assertIn("singleCellSolver", text)


class TestContextResolution(unittest.TestCase):
    """The builder must expose its context-resolution step so tests and
    callers can introspect what slot_keys + values the pipeline will use."""

    def test_resolve_context_collapses_selectors_and_overrides(self) -> None:
        from openfoam_driver.specs.dict_builder import resolve_context
        ctx = resolve_context(
            selectors={"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP"},
            overrides={
                "$ELECTRO_MODEL_COEFFS.solutionAlgorithm": "implicit",
                "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude": "60",
            },
        )
        # Selectors land at their raw key; overrides are stripped of the
        # $ELECTRO_MODEL_COEFFS prefix to match slot_key convention.
        self.assertEqual(ctx["myocardiumSolver"], "monodomainSolver")
        self.assertEqual(ctx["ionicModel"], "TNNP")
        self.assertEqual(ctx["solutionAlgorithm"], "implicit")
        self.assertEqual(ctx["singleCellStimulus.stim_amplitude"], "60")

    def test_resolve_context_overrides_silent_when_none(self) -> None:
        from openfoam_driver.specs.dict_builder import resolve_context
        ctx = resolve_context(
            selectors={"myocardiumSolver": "singleCellSolver"},
            overrides=None,
        )
        self.assertEqual(ctx, {"myocardiumSolver": "singleCellSolver"})


class TestApplicableEntrySelection(unittest.TestCase):
    """Walks IONIC catalog + dict_entries, returns only entries whose
    applicable_when predicate (if any) matches the resolved context."""

    def test_eikonal_context_excludes_ionic_model_entry(self) -> None:
        """ionicModel carries forbidden_when={myocardiumSolver: eikonalSolver}
        — but the *applicable_when*-based exclusion is the tissue entry,
        which applies only to mono/bi/single-cell solvers. Under eikonal,
        the tissue entry must be filtered out by select_applicable_entries."""
        from openfoam_driver.specs.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "eikonalSolver"},
        )
        entries = select_applicable_entries(ctx)
        paths = {e.driver_path for e in entries}
        self.assertNotIn("$ELECTRO_MODEL_COEFFS.tissue", paths)

    def test_monodomain_context_includes_ionic_model_entry(self) -> None:
        from openfoam_driver.specs.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP"},
        )
        entries = select_applicable_entries(ctx)
        paths = {e.driver_path for e in entries}
        self.assertIn("$ELECTRO_MODEL_COEFFS.ionicModel", paths)
        self.assertIn("$ELECTRO_MODEL_COEFFS.tissue", paths)


class TestValuePopulation(unittest.TestCase):
    """Precedence: explicit override > typical_value (when fallback enabled)
    > omit. Returns a dict slot_key -> value for entries that survived
    the applicability filter."""

    def test_override_wins_over_typical_value(self) -> None:
        from openfoam_driver.specs.dict_builder import (
            populate_values,
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
            overrides={
                "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude": "0.4",
            },
        )
        entries = select_applicable_entries(ctx)
        populated = populate_values(entries, ctx, typical_value_fallback=True)
        # Override value wins, not the typical_value="60" from dict_entries.
        self.assertEqual(populated["singleCellStimulus.stim_amplitude"], "0.4")

    def test_typical_value_fills_when_no_override(self) -> None:
        from openfoam_driver.specs.dict_builder import (
            populate_values,
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
        )
        entries = select_applicable_entries(ctx)
        populated = populate_values(entries, ctx, typical_value_fallback=True)
        # singleCellStimulus.stim_amplitude has typical_value="60" in dict_entries.
        self.assertEqual(populated["singleCellStimulus.stim_amplitude"], "60")

    def test_fallback_disabled_omits_typical_value(self) -> None:
        from openfoam_driver.specs.dict_builder import (
            populate_values,
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
        )
        entries = select_applicable_entries(ctx)
        populated = populate_values(entries, ctx, typical_value_fallback=False)
        # No fallback → entry is absent from the populated dict.
        self.assertNotIn("singleCellStimulus.stim_amplitude", populated)

    def test_selector_values_are_present_in_populated_dict(self) -> None:
        """Selectors are part of the context AND many of them correspond to
        DictEntry paths (myocardiumSolver, ionicModel, tissue). Those entries
        must end up in the populated dict using the selector's own value."""
        from openfoam_driver.specs.dict_builder import (
            populate_values,
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
        )
        entries = select_applicable_entries(ctx)
        populated = populate_values(entries, ctx, typical_value_fallback=True)
        self.assertEqual(populated["myocardiumSolver"], "monodomainSolver")
        self.assertEqual(populated["ionicModel"], "TNNP")
        self.assertEqual(populated["tissue"], "epicardialCells")


class TestRequiredCheck(unittest.TestCase):
    """check_required raises ValueError listing every required+applicable
    entry whose slot is missing from the populated dict. Optional entries
    are silently ignored; inapplicable entries are also ignored (filtered
    earlier by select_applicable_entries)."""

    def test_silent_when_all_required_present(self) -> None:
        from openfoam_driver.specs.dict_builder import (
            check_required,
            populate_values,
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
        )
        entries = select_applicable_entries(ctx)
        populated = populate_values(entries, ctx, typical_value_fallback=True)
        # Should not raise — typical_value fallback fills all required leaves.
        check_required(entries, populated, context=ctx)

    def test_raises_listing_missing_required_paths(self) -> None:
        from openfoam_driver.specs.dict_builder import (
            check_required,
            populate_values,
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
        )
        entries = select_applicable_entries(ctx)
        # Fallback OFF — no typical_value fills happen → required-but-no-override
        # entries are missing.
        populated = populate_values(entries, ctx, typical_value_fallback=False)
        with self.assertRaises(ValueError) as ctx_mgr:
            check_required(entries, populated, context=ctx)
        # Error message must enumerate concrete missing paths.
        msg = str(ctx_mgr.exception)
        self.assertIn("singleCellStimulus", msg)

    def test_optional_unset_entries_do_not_raise(self) -> None:
        """A required=False entry that is absent from the populated dict
        is not a violation, even when no typical_value fallback was used."""
        from openfoam_driver.dict_entries import DictEntry
        from openfoam_driver.specs.dict_builder import check_required

        only_optional = [
            DictEntry(
                driver_path="$ELECTRO_MODEL_COEFFS.optionalField",
                description="x",
                source_refs=("ref.C",),
                required=False,
                phases=frozenset({"physics"}),
            ),
        ]
        # populated dict deliberately empty
        check_required(only_optional, {})


class TestValidatorIntegration(unittest.TestCase):
    """build_electro_properties must pass through validate_run before
    returning, so any output an agent receives is validator-clean."""

    def test_build_raises_on_mutex_violation_via_overrides(self) -> None:
        """Setting both stimulusDuration and stimulusDurationList violates the
        structured mutually_exclusive_with constraint — the builder must
        catch it before returning the synthesised text."""
        from openfoam_driver.specs.dict_builder import build_electro_properties
        with self.assertRaises(ValueError) as ctx:
            build_electro_properties(
                selectors={
                    "myocardiumSolver": "monodomainSolver",
                    "ionicModel": "TNNP",
                    "tissue": "epicardialCells",
                },
                overrides={
                    "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDuration": "0.002",
                    "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDurationList": "(0.002 0.001)",
                },
            )
        self.assertIn("mutually exclusive", str(ctx.exception).lower())

    def test_build_raises_on_forbidden_when_violation(self) -> None:
        """ionicModel under eikonalSolver triggers forbidden_when — builder
        must reject this combination."""
        from openfoam_driver.specs.dict_builder import build_electro_properties
        with self.assertRaises(ValueError) as ctx:
            build_electro_properties(
                selectors={
                    "myocardiumSolver": "eikonalSolver",
                    "ionicModel": "TNNP",
                },
            )
        self.assertIn("forbidden", str(ctx.exception).lower())


class TestSerialisation(unittest.TestCase):
    """The output must be a real OpenFOAM dict: FoamFile preamble, top-level
    selectors at the root, everything else nested under <solver>Coeffs with
    sub-blocks emitted recursively."""

    def test_singlecell_output_has_solver_coeffs_block(self) -> None:
        from openfoam_driver.specs.dict_builder import build_electro_properties
        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "singleCellSolver",
                "ionicModel": "AlievPanfilov",
                "tissue": "myocyte",
            },
        )
        self.assertIn("myocardiumSolver singleCellSolver;", text)
        self.assertIn("singleCellSolverCoeffs", text)
        self.assertIn("ionicModel AlievPanfilov;", text)
        self.assertIn("tissue myocyte;", text)

    def test_singlecell_output_nests_stimulus_subblock(self) -> None:
        from openfoam_driver.specs.dict_builder import build_electro_properties
        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "singleCellSolver",
                "ionicModel": "AlievPanfilov",
                "tissue": "myocyte",
            },
        )
        # singleCellStimulus is a sub-block with its own children
        self.assertIn("singleCellStimulus", text)
        # The stim_amplitude leaf must appear in nested form (typical_value="60")
        self.assertIn("stim_amplitude 60;", text)
        self.assertIn("stim_start 0.0;", text)

    def test_monodomain_output_uses_monodomain_coeffs(self) -> None:
        from openfoam_driver.specs.dict_builder import build_electro_properties
        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "monodomainSolver",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        )
        self.assertIn("myocardiumSolver monodomainSolver;", text)
        self.assertIn("monodomainSolverCoeffs", text)
        self.assertNotIn("singleCellSolverCoeffs", text)

    def test_output_parses_back_with_existing_helpers(self) -> None:
        """Round-trip: build → write to disk → re-parse with detect_* helpers
        → confirm input matches what's recovered."""
        import tempfile
        from pathlib import Path
        from openfoam_driver.specs.dict_builder import build_electro_properties
        from openfoam_driver.specs.common import (
            detect_ionic_model_name,
            detect_myocardium_solver_name,
        )
        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "monodomainSolver",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        )
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "electroProperties"
            path.write_text(text)
            self.assertEqual(detect_myocardium_solver_name(path), "monodomainSolver")
            self.assertEqual(detect_ionic_model_name(path), "TNNP")


class TestPhysicsPropertiesBuilder(unittest.TestCase):
    """build_physics_properties mirrors the electroProperties pipeline
    against the small PHYSICS_PROPERTY_ENTRIES set. No <solver>Coeffs
    wrapper — physics keys live at the dict root."""

    def test_function_accepts_documented_kwargs(self) -> None:
        import inspect
        from openfoam_driver.specs.dict_builder import build_physics_properties
        sig = inspect.signature(build_physics_properties)
        params = sig.parameters
        self.assertIn("selectors", params)
        self.assertIn("overrides", params)
        self.assertIn("typical_value_fallback", params)
        self.assertEqual(params["overrides"].kind, inspect.Parameter.KEYWORD_ONLY)

    def test_minimal_electroModel_build(self) -> None:
        from openfoam_driver.specs.dict_builder import build_physics_properties
        text = build_physics_properties(selectors={"type": "electroModel"})
        self.assertIn("type electroModel;", text)
        self.assertIn("FoamFile", text)
        # Preamble must point at physicsProperties, not electroProperties.
        self.assertIn("object      physicsProperties", text)
        self.assertNotIn("electroProperties", text)
        # No <solver>Coeffs block — physics keys are root-level.
        self.assertNotIn("Coeffs", text)

    def test_missing_required_type_raises(self) -> None:
        from openfoam_driver.specs.dict_builder import build_physics_properties
        with self.assertRaises(ValueError) as ctx:
            build_physics_properties(selectors={})
        self.assertIn("type", str(ctx.exception))

    def test_invalid_enum_value_raises(self) -> None:
        from openfoam_driver.specs.dict_builder import build_physics_properties
        with self.assertRaises(ValueError) as ctx:
            build_physics_properties(selectors={"type": "notARealModel"})
        self.assertIn("notARealModel", str(ctx.exception))


class TestBuildAndLaunch(unittest.TestCase):
    """build_and_launch closes the last gap between "agent can construct
    a dict" and "agent can launch a run". It writes both dicts to a case
    directory and invokes the engine via the generic_case spec factory.
    """

    def test_writes_both_dicts_to_case_dir(self) -> None:
        """The dry_run=True path writes the dicts and exits without
        running cardiacFoam, so the test never needs the binary."""
        import tempfile
        from pathlib import Path
        from openfoam_driver.specs.dict_builder import build_and_launch

        with tempfile.TemporaryDirectory() as temp:
            case_dir = Path(temp) / "case"
            result = build_and_launch(
                electro_selectors={
                    "myocardiumSolver": "singleCellSolver",
                    "ionicModel": "AlievPanfilov",
                    "tissue": "myocyte",
                },
                physics_selectors={"type": "electroModel"},
                case_dir=case_dir,
                dry_run=True,
            )
            self.assertTrue((case_dir / "constant" / "electroProperties").exists())
            self.assertTrue((case_dir / "constant" / "physicsProperties").exists())
            self.assertEqual(result["case_dir"], str(case_dir))
            self.assertEqual(result["status"], "dry_run_complete")

    def test_existing_case_dir_is_not_overwritten_without_consent(self) -> None:
        """The wrapper must refuse to clobber an existing case_dir unless
        the caller explicitly passes `overwrite=True`."""
        import tempfile
        from pathlib import Path
        from openfoam_driver.specs.dict_builder import build_and_launch

        with tempfile.TemporaryDirectory() as temp:
            case_dir = Path(temp) / "case"
            (case_dir / "constant").mkdir(parents=True)
            (case_dir / "constant" / "electroProperties").write_text("# pre-existing\n")

            with self.assertRaises(FileExistsError):
                build_and_launch(
                    electro_selectors={
                        "myocardiumSolver": "singleCellSolver",
                        "ionicModel": "AlievPanfilov",
                        "tissue": "myocyte",
                    },
                    physics_selectors={"type": "electroModel"},
                    case_dir=case_dir,
                    dry_run=True,
                )

    def test_overwrite_true_replaces_existing_dicts(self) -> None:
        import tempfile
        from pathlib import Path
        from openfoam_driver.specs.dict_builder import build_and_launch

        with tempfile.TemporaryDirectory() as temp:
            case_dir = Path(temp) / "case"
            (case_dir / "constant").mkdir(parents=True)
            old_text = "# pre-existing electroProperties\n"
            (case_dir / "constant" / "electroProperties").write_text(old_text)

            build_and_launch(
                electro_selectors={
                    "myocardiumSolver": "singleCellSolver",
                    "ionicModel": "AlievPanfilov",
                    "tissue": "myocyte",
                },
                physics_selectors={"type": "electroModel"},
                case_dir=case_dir,
                dry_run=True,
                overwrite=True,
            )
            text = (case_dir / "constant" / "electroProperties").read_text()
            self.assertNotEqual(text, old_text)
            self.assertIn("myocardiumSolver singleCellSolver;", text)


class TestBuildAndLaunchDirectRun(unittest.TestCase):
    """build_and_launch passes solver_command='cardiacFoam' to make_spec."""

    def _base_selectors(self):
        return (
            {"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
            {"type": "electroModel"},
        )

    def test_dry_run_still_completes_without_pre_solve(self) -> None:
        import tempfile
        from pathlib import Path
        from openfoam_driver.specs.dict_builder import build_and_launch
        electro, physics = self._base_selectors()
        with tempfile.TemporaryDirectory() as d:
            result = build_and_launch(
                electro,
                physics_selectors=physics,
                case_dir=Path(d) / "case",
                dry_run=True,
                pre_solve_commands=["vtkUnstructuredToFoam"],
            )
        self.assertEqual(result["status"], "dry_run_complete")

    def test_pre_solve_commands_run_before_solver(self) -> None:
        import subprocess
        import tempfile
        from pathlib import Path
        from unittest.mock import patch
        from openfoam_driver.specs.dict_builder import build_and_launch
        electro, physics = self._base_selectors()
        with tempfile.TemporaryDirectory() as d:
            case_dir = Path(d) / "case"
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = subprocess.CompletedProcess([], 0)
                build_and_launch(
                    electro,
                    physics_selectors=physics,
                    case_dir=case_dir,
                    pre_solve_commands=["vtkUnstructuredToFoam"],
                )
            calls = mock_run.call_args_list
            self.assertGreaterEqual(len(calls), 2)
            first_args = calls[0].args[0]
            solver_args = calls[-1].args[0]
            self.assertIn("vtkUnstructuredToFoam", first_args)
            self.assertIn("cardiacFoam", solver_args)

    def test_no_pre_solve_calls_only_solver(self) -> None:
        import subprocess
        import tempfile
        from pathlib import Path
        from unittest.mock import patch
        from openfoam_driver.specs.dict_builder import build_and_launch
        electro, physics = self._base_selectors()
        with tempfile.TemporaryDirectory() as d:
            case_dir = Path(d) / "case"
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = subprocess.CompletedProcess([], 0)
                build_and_launch(
                    electro,
                    physics_selectors=physics,
                    case_dir=case_dir,
                )
            calls = mock_run.call_args_list
            self.assertEqual(len(calls), 1)
            self.assertIn("cardiacFoam", calls[0].args[0])


if __name__ == "__main__":
    unittest.main()
