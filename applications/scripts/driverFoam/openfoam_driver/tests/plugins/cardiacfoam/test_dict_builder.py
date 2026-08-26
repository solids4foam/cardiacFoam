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
#     test_dict_builder
#
# Description
#     Tests dict builder logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the dict_builder.

The builder synthesizes a complete electroProperties dict from
minimum-viable agent intent (selectors + overrides). It enforces the same
constraints the validator does, so its output is guaranteed validator-clean
by construction.
"""
from __future__ import annotations

import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[7]
SINGLE_CELL_ELECTRO_PROPERTIES = (
    REPO_ROOT / "tutorials" / "electrophysiologyProtocols" / "singleCell"
    / "constant" / "electroProperties"
)
PURKINJE_NIEDERER = REPO_ROOT / "tutorials" / "NiedererEtAl2011" / "purkinjeNiedererEtAl2011"
PURKINJE_ELECTRO_PROPERTIES_MONODOMAIN = (
    PURKINJE_NIEDERER / "constant" / "electroProperties.monodomain"
)


class TestDictBuilderModule(unittest.TestCase):
    """Module-level structural contract — the import path + signature."""

    def test_module_exposes_build_electro_properties(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
        self.assertTrue(callable(build_electro_properties))

    def test_function_accepts_documented_kwargs(self) -> None:
        import inspect
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import resolve_context
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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import resolve_context
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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
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

    def test_field_source_omits_monodomain_uniform_tensor(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties

        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "monodomainSolver",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
                "conductivitySource": "field",
            },
        )
        self.assertIn("conductivitySource field;", text)
        self.assertNotIn("\n    conductivity ", text)

    def test_field_source_omits_both_bidomain_uniform_tensors(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties

        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "bidomainSolver",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
                "conductivitySource": "field",
            },
        )
        self.assertIn("conductivitySource field;", text)
        self.assertNotIn("conductivityIntracellular", text)
        self.assertNotIn("conductivityExtracellular", text)

    def test_spatial_solver_defaults_to_uniform_source(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties

        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "monodomainSolver",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        )
        self.assertIn("conductivitySource uniform;", text)
        self.assertIn("\n    conductivity ", text)


class TestValuePopulation(unittest.TestCase):
    """Precedence: explicit override > typical_value (when fallback enabled)
    > omit. Returns a dict slot_key -> value for entries that survived
    the applicability filter."""

    def test_override_wins_over_typical_value(self) -> None:
        from openfoam_driver.core.specs.dict_builder import populate_values
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
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
        from openfoam_driver.core.specs.dict_builder import populate_values
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )
        # A stimulus override makes the block "configured", which is what the
        # family is gated on -- an unconfigured case must NOT get a stimulus
        # invented for it (stimulusIO.C:149-155 treats an absent block as a
        # legal no-op protocol).
        ctx = resolve_context(
            selectors={"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
            overrides={"$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_start": "20"},
        )
        entries = select_applicable_entries(ctx)
        populated = populate_values(entries, ctx, typical_value_fallback=True)
        # singleCellStimulus.stim_amplitude has typical_value="60" in dict_entries.
        self.assertEqual(populated["singleCellStimulus.stim_amplitude"], "60")

    def test_fallback_disabled_omits_typical_value(self) -> None:
        from openfoam_driver.core.specs.dict_builder import populate_values
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
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
        from openfoam_driver.core.specs.dict_builder import populate_values
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
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
        from openfoam_driver.core.specs.dict_builder import (
            check_required,
            populate_values,
        )
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
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
        from openfoam_driver.core.specs.dict_builder import (
            check_required,
            populate_values,
        )
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )
        ctx = resolve_context(
            selectors={"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
            overrides={"$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_start": "20"},
        )
        entries = select_applicable_entries(ctx)
        # Fallback OFF — no typical_value fills happen → required-but-no-override
        # entries are missing. The stimulus override above makes the stimulus
        # family applicable, so its four guarded keys appear in the listing.
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
        from openfoam_driver.core.specs.dict_builder import check_required

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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
        with self.assertRaises(ValueError) as ctx:
            build_electro_properties(
                selectors={
                    "myocardiumSolver": "eikonalSolver",
                    "ionicModel": "TNNP",
                },
            )
        self.assertIn("forbidden", str(ctx.exception).lower())


class TestSerialisation(unittest.TestCase):
    """The output must be a real OpenFOAM dict. We use a snapshot test
    to verify the generated electroProperties output matches expectations."""

    def test_singlecell_output_matches_snapshot(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
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
        # No stimulus was asked for, so none is invented. stimulusIO.C:149-155
        # treats an absent singleCellStimulus block as a legal no-op protocol;
        # filling it from typical_value would silently pace a quiescent case.
        self.assertNotIn("singleCellStimulus", text)
        self.assertNotIn("stim_amplitude", text)

    def test_singlecell_stimulus_appears_once_configured(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
        text = build_electro_properties(
            selectors={
                "myocardiumSolver": "singleCellSolver",
                "ionicModel": "AlievPanfilov",
                "tissue": "myocyte",
            },
            overrides={
                "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_start": "20",
            },
        )
        # One override configures the block; the rest of the family then fills
        # from typical_value, so the four keys stimulusIO.C:159-176 requires
        # together are never half-written.
        self.assertIn("singleCellStimulus", text)
        self.assertIn("stim_start 20;", text)
        self.assertIn("stim_amplitude 60;", text)

    def test_monodomain_output_matches_snapshot(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
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


class TestPhysicsPropertiesBuilder(unittest.TestCase):
    """build_physics_properties mirrors the electroProperties pipeline
    against the small PHYSICS_PROPERTY_ENTRIES set. No <solver>Coeffs
    wrapper — physics keys live at the dict root."""

    def test_function_accepts_documented_kwargs(self) -> None:
        import inspect
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_physics_properties
        sig = inspect.signature(build_physics_properties)
        params = sig.parameters
        self.assertIn("selectors", params)
        self.assertIn("overrides", params)
        self.assertIn("typical_value_fallback", params)
        self.assertEqual(params["overrides"].kind, inspect.Parameter.KEYWORD_ONLY)

    def test_minimal_electroModel_build(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_physics_properties
        text = build_physics_properties(selectors={"type": "electroModel"})
        self.assertIn("type electroModel;", text)
        self.assertIn("FoamFile", text)
        # Preamble must point at physicsProperties, not electroProperties.
        self.assertIn("object      physicsProperties", text)
        self.assertNotIn("electroProperties", text)
        # No <solver>Coeffs block — physics keys are root-level.
        self.assertNotIn("Coeffs", text)

    def test_missing_required_type_raises(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_physics_properties
        with self.assertRaises(ValueError) as ctx:
            build_physics_properties(selectors={})
        self.assertIn("type", str(ctx.exception))

    def test_invalid_enum_value_raises(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_physics_properties
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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch

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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch

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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch

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


class TestBuildAndLaunchMeshProvisioning(unittest.TestCase):
    """build_and_launch must leave every case with a real mesh on disk.

    electroModel.C requires a real fvMesh regardless of solver (confirmed via
    `refCast<const fvMesh>(mesh())` at electroModel.C:344) -- even
    singleCellSolver needs one. Neither build_and_launch nor
    sweep_runner.materialize_case provisioned any mesh before this fix, so no
    solver built from scratch via sweep-run/case_folder could ever complete
    (see project_driverfoam_sweep_bugs_found memory item #3).
    """

    def test_single_cell_solver_gets_a_static_polymesh(self) -> None:
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch

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
            poly_mesh = case_dir / "constant" / "polyMesh"
            for name in ("points", "faces", "owner", "neighbour", "boundary"):
                self.assertTrue((poly_mesh / name).exists(), f"missing {name}")
            self.assertFalse(result.get("needs_block_mesh", False))

    def test_spatial_solver_gets_a_block_mesh_dict(self) -> None:
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch

        with tempfile.TemporaryDirectory() as temp:
            case_dir = Path(temp) / "case"
            result = build_and_launch(
                electro_selectors={
                    "myocardiumSolver": "monodomainSolver",
                    "ionicModel": "TNNP",
                    "tissue": "epicardialCells",
                },
                physics_selectors={"type": "electroModel"},
                case_dir=case_dir,
                dry_run=True,
            )
            block_mesh_dict = case_dir / "system" / "blockMeshDict"
            self.assertTrue(block_mesh_dict.exists())
            self.assertIn("blocks", block_mesh_dict.read_text())
            self.assertFalse((case_dir / "constant" / "polyMesh").exists())
            self.assertTrue(result.get("needs_block_mesh", False))

    def test_dx_kwarg_controls_generated_block_mesh_resolution(self) -> None:
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch
        from openfoam_driver.core.specs.mesh_provisioning import default_block_mesh_dict_text

        with tempfile.TemporaryDirectory() as temp:
            case_dir = Path(temp) / "case"
            build_and_launch(
                electro_selectors={
                    "myocardiumSolver": "monodomainSolver",
                    "ionicModel": "TNNP",
                    "tissue": "epicardialCells",
                },
                physics_selectors={"type": "electroModel"},
                case_dir=case_dir,
                dry_run=True,
                dx=0.0004,
            )
            written = (case_dir / "system" / "blockMeshDict").read_text()
            self.assertEqual(written, default_block_mesh_dict_text(dx_m=0.0004))
            self.assertNotEqual(written, default_block_mesh_dict_text())

    def test_dx_kwarg_rejected_for_meshless_solver(self) -> None:
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch

        with tempfile.TemporaryDirectory() as temp:
            case_dir = Path(temp) / "case"
            with self.assertRaisesRegex(ValueError, "dx"):
                build_and_launch(
                    electro_selectors={
                        "myocardiumSolver": "singleCellSolver",
                        "ionicModel": "AlievPanfilov",
                        "tissue": "myocyte",
                    },
                    physics_selectors={"type": "electroModel"},
                    case_dir=case_dir,
                    dry_run=True,
                    dx=0.0004,
                )

    def test_existing_mesh_is_not_clobbered_without_overwrite(self) -> None:
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch

        with tempfile.TemporaryDirectory() as temp:
            case_dir = Path(temp) / "case"
            block_mesh_dict = case_dir / "system" / "blockMeshDict"
            block_mesh_dict.parent.mkdir(parents=True)
            block_mesh_dict.write_text("// pre-existing custom mesh\n")
            (case_dir / "constant").mkdir(parents=True)
            (case_dir / "constant" / "electroProperties").write_text("# pre-existing\n")
            build_and_launch(
                electro_selectors={
                    "myocardiumSolver": "monodomainSolver",
                    "ionicModel": "TNNP",
                    "tissue": "epicardialCells",
                },
                physics_selectors={"type": "electroModel"},
                case_dir=case_dir,
                dry_run=True,
                overwrite=True,
            )
            self.assertEqual(block_mesh_dict.read_text(), "// pre-existing custom mesh\n")


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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch
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
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch
        electro, physics = self._base_selectors()
        with tempfile.TemporaryDirectory() as d:
            case_dir = Path(d) / "case"
            with patch("subprocess.run") as mock_run:
                # capture_output=True, text=True always yields str stdout/stderr;
                # CompletedProcess defaults them to None, which is not a
                # shape subprocess.run can actually return. The global patch
                # is intentional (the test observes calls from several
                # modules), so the mock must be faithful instead.
                mock_run.return_value = subprocess.CompletedProcess(
                    [], 0, stdout="", stderr=""
                )
                build_and_launch(
                    electro,
                    physics_selectors=physics,
                    case_dir=case_dir,
                    pre_solve_commands=["vtkUnstructuredToFoam"],
                )
            calls = mock_run.call_args_list
            # calls[0] is now the strict path's own environment load (sources
            # the OpenFOAM bashrc to capture env vars, same as every other
            # run --strict invocation does) -- find the actual step calls by
            # content rather than a fixed index.
            self.assertGreaterEqual(len(calls), 2)
            pre_solve_index = next(i for i, c in enumerate(calls) if "vtkUnstructuredToFoam" in c.args[0])
            solver_index = next(i for i, c in enumerate(calls) if "cardiacFoam" in c.args[0])
            self.assertLess(pre_solve_index, solver_index)

    def test_no_pre_solve_calls_only_solver(self) -> None:
        import subprocess
        import tempfile
        from pathlib import Path
        from unittest.mock import patch
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch
        electro, physics = self._base_selectors()
        with tempfile.TemporaryDirectory() as d:
            case_dir = Path(d) / "case"
            with patch("subprocess.run") as mock_run:
                # capture_output=True, text=True always yields str stdout/stderr;
                # CompletedProcess defaults them to None, which is not a
                # shape subprocess.run can actually return. The global patch
                # is intentional (the test observes calls from several
                # modules), so the mock must be faithful instead.
                mock_run.return_value = subprocess.CompletedProcess(
                    [], 0, stdout="", stderr=""
                )
                build_and_launch(
                    electro,
                    physics_selectors=physics,
                    case_dir=case_dir,
                )
            calls = mock_run.call_args_list
            self.assertGreaterEqual(len(calls), 1)
            self.assertIn("cardiacFoam", calls[-1].args[0])


class TestParseElectroProperties(unittest.TestCase):
    """parse_electro_properties reads an existing file back into selectors+overrides."""

    @staticmethod
    def _build_and_write(tmp_dir, selectors, overrides=None):
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
        text = build_electro_properties(selectors, overrides=overrides)
        p = Path(tmp_dir) / "electroProperties"
        p.write_text(text)
        return p

    def test_returns_dict_with_selectors_and_overrides_keys(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
        with tempfile.TemporaryDirectory() as d:
            p = self._build_and_write(
                d,
                {"myocardiumSolver": "singleCellSolver",
                 "ionicModel": "AlievPanfilov",
                 "tissue": "myocyte"},
            )
            result = parse_electro_properties(p)
            self.assertIn("selectors", result)
            self.assertIn("overrides", result)

    def test_solver_in_selectors(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
        with tempfile.TemporaryDirectory() as d:
            p = self._build_and_write(
                d,
                {"myocardiumSolver": "monodomainSolver",
                 "ionicModel": "TNNP",
                 "tissue": "epicardialCells"},
            )
            result = parse_electro_properties(p)
            self.assertEqual(result["selectors"]["myocardiumSolver"], "monodomainSolver")

    def test_ionic_model_and_tissue_in_selectors(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
        with tempfile.TemporaryDirectory() as d:
            p = self._build_and_write(
                d,
                {"myocardiumSolver": "singleCellSolver",
                 "ionicModel": "AlievPanfilov",
                 "tissue": "myocyte"},
            )
            result = parse_electro_properties(p)
            self.assertEqual(result["selectors"]["ionicModel"], "AlievPanfilov")
            self.assertEqual(result["selectors"]["tissue"], "myocyte")

    def test_ignored_keys_lists_structurally_skipped_dynamic_paths(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
        with tempfile.TemporaryDirectory() as d:
            p = self._build_and_write(
                d,
                {"myocardiumSolver": "monodomainSolver",
                 "ionicModel": "TNNP",
                 "tissue": "epicardialCells"},
            )
            result = parse_electro_properties(p)
            # Backward compatible: selectors/overrides still present.
            self.assertIn("selectors", result)
            self.assertIn("overrides", result)
            # New: the parser now surfaces the driver_path families it does not
            # round-trip instead of dropping them silently.
            self.assertIn("ignored_keys", result)
            self.assertIsInstance(result["ignored_keys"], list)
            # dynamic_path entries exist in the catalog, so the list is non-empty
            # and every entry is a $ELECTRO_MODEL_COEFFS driver_path.
            self.assertTrue(result["ignored_keys"])
            self.assertTrue(
                all(k.startswith("$ELECTRO_MODEL_COEFFS") for k in result["ignored_keys"])
            )

    def test_non_default_override_is_captured(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
        with tempfile.TemporaryDirectory() as d:
            p = self._build_and_write(
                d,
                {"myocardiumSolver": "monodomainSolver",
                 "ionicModel": "AlievPanfilov",
                 "tissue": "myocyte"},
                overrides={
                    "$ELECTRO_MODEL_COEFFS.solutionAlgorithm": "implicit",
                },
            )
            result = parse_electro_properties(p)
            self.assertEqual(
                result["overrides"].get("$ELECTRO_MODEL_COEFFS.solutionAlgorithm"),
                "implicit",
            )

    def test_default_value_absent_from_overrides(self) -> None:
        """typical_value entries that were not changed must not appear as overrides."""
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
        with tempfile.TemporaryDirectory() as d:
            p = self._build_and_write(
                d,
                {"myocardiumSolver": "singleCellSolver",
                 "ionicModel": "AlievPanfilov",
                 "tissue": "myocyte"},
            )
            result = parse_electro_properties(p)
            self.assertNotIn(
                "$ELECTRO_MODEL_COEFFS.solutionAlgorithm",
                result["overrides"],
            )

    def test_selector_keys_not_duplicated_in_overrides(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
        with tempfile.TemporaryDirectory() as d:
            p = self._build_and_write(
                d,
                {"myocardiumSolver": "singleCellSolver",
                 "ionicModel": "AlievPanfilov",
                 "tissue": "myocyte"},
            )
            result = parse_electro_properties(p)
            override_slot_keys = {
                k.replace("$ELECTRO_MODEL_COEFFS.", "")
                for k in result["overrides"]
            }
            for sel_key in ("myocardiumSolver", "ionicModel", "tissue"):
                self.assertNotIn(sel_key, override_slot_keys)

    def test_active_tension_model_survives_singlecell_roundtrip(self) -> None:
        """Regression test: activeTensionModel is a flat word entry directly
        inside singleCellSolverCoeffs (see singleCellSolver.C's
        electroProperties().found("activeTensionModel")) — not a
        'activeTensionModel { activeTensionModel <x>; }' sub-block. Both
        build_electro_properties (synthesis) and parse_electro_properties
        (round-trip) must preserve it for singleCellSolver."""
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            parse_electro_properties,
        )
        selectors = {
            "myocardiumSolver": "singleCellSolver",
            "ionicModel": "TWorld",
            "tissue": "endocardialCells",
        }
        overrides = {"$ELECTRO_MODEL_COEFFS.activeTensionModel": "LandNiederer"}

        text = build_electro_properties(selectors, overrides=overrides)
        self.assertIn("activeTensionModel LandNiederer;", text)
        self.assertNotIn("activeTensionModel\n    {", text)

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "electroProperties"
            p.write_text(text)
            parsed = parse_electro_properties(p)

        self.assertEqual(
            parsed["overrides"].get("$ELECTRO_MODEL_COEFFS.activeTensionModel"),
            "LandNiederer",
        )

    def test_nested_dynamic_placeholders_survive_roundtrip(self) -> None:
        """Expand both the tissue scope and constant-name placeholders.

        A one-placeholder parser can discover ``global`` but then looks for
        the literal ``<constant_name>`` key, silently dropping the override.
        """
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            parse_electro_properties,
        )

        driver_path = "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.set.R"
        text = build_electro_properties(
            {
                "myocardiumSolver": "monodomainSolver",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
            overrides={driver_path: "1.23"},
        )
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "electroProperties"
            path.write_text(text)
            parsed = parse_electro_properties(path)

        self.assertEqual(parsed["overrides"].get(driver_path), "1.23")
        self.assertNotIn(
            "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.<scope>.set.<constant_name>",
            parsed["ignored_keys"],
        )

    def test_active_tension_model_recovered_from_real_singlecell_tutorial(self) -> None:
        """The hand-authored singleCell tutorial dict declares
        'activeTensionModel LandNiederer;' as a flat entry — parsing it must
        not silently drop that setting."""
        from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties

        if not SINGLE_CELL_ELECTRO_PROPERTIES.exists():
            self.skipTest("tutorial fixture not present in this checkout")

        result = parse_electro_properties(SINGLE_CELL_ELECTRO_PROPERTIES)
        self.assertEqual(
            result["overrides"].get("$ELECTRO_MODEL_COEFFS.activeTensionModel"),
            "LandNiederer",
        )

    def test_roundtrip_produces_equivalent_text(self) -> None:
        """build → write → parse → rebuild must produce a semantically
        equivalent FoamFile: re-parsing the rebuilt text must yield the same
        selectors/overrides as the original. We deliberately don't compare
        raw text here — `parse_electro_properties` reads values back via
        `foamDictionary`, which canonicalizes numeric formatting as a side
        effect of parsing (e.g. `0.0` -> `0`, `1e-6` -> `1e-06`), so the
        rebuilt text can legitimately differ cosmetically from the original
        while still describing the same dict.
        """
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            parse_electro_properties,
        )
        original_selectors = {
            "myocardiumSolver": "monodomainSolver",
            "ionicModel": "TNNP",
            "tissue": "epicardialCells",
        }
        original_overrides = {
            "$ELECTRO_MODEL_COEFFS.solutionAlgorithm": "explicit",
        }
        original_text = build_electro_properties(
            original_selectors, overrides=original_overrides,
        )
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "electroProperties"
            p.write_text(original_text)
            parsed = parse_electro_properties(p)

        rebuilt_text = build_electro_properties(
            parsed["selectors"], overrides=parsed["overrides"] or None,
        )
        with tempfile.TemporaryDirectory() as d2:
            p2 = Path(d2) / "electroProperties"
            p2.write_text(rebuilt_text)
            reparsed = parse_electro_properties(p2)

        self.assertEqual(parsed, reparsed)


class TestBuildAndLaunchControlDict(unittest.TestCase):
    """build_and_launch delta_t / end_time patch an existing system/controlDict."""

    @staticmethod
    def _selectors():
        return (
            {"myocardiumSolver": "singleCellSolver",
             "ionicModel": "AlievPanfilov",
             "tissue": "myocyte"},
            {"type": "electroModel"},
        )

    @staticmethod
    def _make_case_with_control_dict(d: str) -> "object":
        from pathlib import Path
        case_dir = Path(d) / "case"
        (case_dir / "system").mkdir(parents=True)
        (case_dir / "system" / "controlDict").write_text(
            "deltaT    0.05;\nendTime   1.0;\n"
        )
        return case_dir

    def test_delta_t_written_to_control_dict(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch
        electro, physics = self._selectors()
        with tempfile.TemporaryDirectory() as d:
            case_dir = self._make_case_with_control_dict(d)
            build_and_launch(
                electro,
                physics_selectors=physics,
                case_dir=case_dir,
                delta_t=0.001,
                dry_run=True,
            )
            text = (case_dir / "system" / "controlDict").read_text()
            self.assertIn("0.001", text)

    def test_end_time_written_to_control_dict(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch
        electro, physics = self._selectors()
        with tempfile.TemporaryDirectory() as d:
            case_dir = self._make_case_with_control_dict(d)
            build_and_launch(
                electro,
                physics_selectors=physics,
                case_dir=case_dir,
                end_time=0.002,
                dry_run=True,
            )
            text = (case_dir / "system" / "controlDict").read_text()
            self.assertIn("0.002", text)

    def test_none_params_leave_control_dict_unchanged(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch
        electro, physics = self._selectors()
        with tempfile.TemporaryDirectory() as d:
            case_dir = self._make_case_with_control_dict(d)
            original = (case_dir / "system" / "controlDict").read_text()
            build_and_launch(
                electro,
                physics_selectors=physics,
                case_dir=case_dir,
                dry_run=True,
            )
            self.assertEqual(
                (case_dir / "system" / "controlDict").read_text(),
                original,
            )

    def test_control_dict_is_generated_when_delta_t_set(self) -> None:
        import tempfile
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_and_launch
        electro, physics = self._selectors()
        with tempfile.TemporaryDirectory() as d:
            case_dir = Path(d) / "case"
            case_dir.mkdir()
            build_and_launch(
                electro,
                physics_selectors=physics,
                case_dir=case_dir,
                delta_t=0.001,
                dry_run=True,
            )
            control_dict_path = case_dir / "system" / "controlDict"
            self.assertTrue(control_dict_path.exists())
            self.assertIn("deltaT", control_dict_path.read_text())


class TestEikonalECGHeterogeneity(unittest.TestCase):
    """sigmaExtracellular and ionicHeterogeneity entries for eikonalSolver +
    eikonalECG: catalog visibility, applicable_when firing, and round-trip."""

    def test_sigmaExtracellular_in_catalog_when_ecgDomains_present(self) -> None:
        """sigmaExtracellular must appear in select_applicable_entries whenever
        any ecgDomains override is set ($ecgDomains_present virtual key)."""
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )
        from openfoam_driver.core.specs.validation import slot_key

        context = resolve_context(
            selectors={"myocardiumSolver": "eikonalSolver"},
            overrides={
                "$ELECTRO_MODEL_COEFFS.ecgDomains.ECG.ecgSolver": "eikonalECG",
            },
        )
        entries = select_applicable_entries(context)
        paths = {e.driver_path for e in entries}
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.sigmaExtracellular",
            paths,
        )

    def test_sigmaExtracellular_absent_without_ecgDomains(self) -> None:
        """sigmaExtracellular must NOT appear when no ecgDomains are configured."""
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )

        context = resolve_context(selectors={"myocardiumSolver": "eikonalSolver"})
        entries = select_applicable_entries(context)
        paths = {e.driver_path for e in entries}
        self.assertNotIn(
            "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.sigmaExtracellular",
            paths,
        )

    def test_ionic_heterogeneity_entries_applicable_for_eikonalSolver(self) -> None:
        """All ionicHeterogeneity sub-entries must be selectable for eikonalSolver."""
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )

        context = resolve_context(selectors={"myocardiumSolver": "eikonalSolver"})
        entries = select_applicable_entries(context)
        paths = {e.driver_path for e in entries}
        expected = {
            "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.field",
            "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mode",
            "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.endoMInterface",
            "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mEpiInterface",
            "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionWidth",
            "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionMode",
            "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.smoothing",
        }
        self.assertTrue(expected.issubset(paths), f"Missing: {expected - paths}")

    def test_ionic_heterogeneity_entries_applicable_for_monodomainSolver(self) -> None:
        """ionicHeterogeneity entries must still fire for monodomainSolver (no regression)."""
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )

        context = resolve_context(
            selectors={"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
        )
        entries = select_applicable_entries(context)
        paths = {e.driver_path for e in entries}
        self.assertIn("$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionMode", paths)
        self.assertIn("$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.smoothing", paths)

    def test_eikonalSolver_with_ionicHeterogeneity_overrides_round_trips(self) -> None:
        """build_electro_properties must accept eikonalSolver with an explicit
        ionicHeterogeneity block and write the values into the output dict."""
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties

        text = build_electro_properties(
            selectors={"myocardiumSolver": "eikonalSolver"},
            overrides={
                # Minimum required eikonalSolver fields without typical_values
                "$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach": "true",
                "$ELECTRO_MODEL_COEFFS.stimulusLocationMin": "(0 0 0)",
                "$ELECTRO_MODEL_COEFFS.stimulusLocationMax": "(0.01 0.01 0.01)",
                "$ELECTRO_MODEL_COEFFS.c0": "60",
                # ionicHeterogeneity block for eikonalECG blend mode
                "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.field": "uvc_transmural",
                "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mode": "transmuralBands",
                "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.endoMInterface": "0.3",
                "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mEpiInterface": "0.7",
                "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionWidth": "0.1",
                "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionMode": "blend",
                "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.smoothing": "smoothstep",
            },
        )
        self.assertIn("ionicHeterogeneity", text)
        self.assertIn("transitionMode blend;", text)
        self.assertIn("transitionWidth 0.1;", text)
        self.assertIn("smoothing smoothstep;", text)
        self.assertIn("field uvc_transmural;", text)


class TestRestitutionEikonalSolver1D(unittest.TestCase):
    """Regression coverage for the conductionSystemSolver applicable_when
    bug: selecting restitutionEikonalSolver1D must emit its own
    solver-specific keys (useEdgeConductance, referenceConductance), not a
    dict that is byte-identical to plain eikonalSolver1D apart from the
    solver name.

    Root cause was that select_applicable_entries never resolved the
    "$SCOPE." prefix or the "<name>" placeholder on applicable_when
    predicate keys before comparing them against the (prefix-stripped,
    instance-resolved) context -- so any applicable_when that referenced a
    real driver_path instead of a bare virtual "$..._present" token could
    never match, silently dropping the gated entry regardless of context.
    See openfoam_driver/specs/validation.py's _predicate_matches.
    """

    _NETWORK_NAME = "purkinjeNetwork"
    _COUPLING_NAME = "pvj"

    @classmethod
    def _overrides(cls, conduction_system_solver: str) -> dict[str, str]:
        purkinje = (
            f"$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.{cls._NETWORK_NAME}"
            ".purkinjeGraphModelCoeffs"
        )
        coupling = (
            f"$ELECTRO_MODEL_COEFFS.domainCouplings.{cls._COUPLING_NAME}"
        )
        return {
            f"$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.{cls._NETWORK_NAME}"
            ".conductionSystemDomain": "purkinjeGraphModel",
            f"{purkinje}.graphFile": "purkinjeGraph",
            f"{purkinje}.purkinjeCV": "[0 1 -1 0 0 0 0] 2.0",
            f"{purkinje}.vm1DRest": "-0.084",
            f"{purkinje}.rootStimulus.node": "0",
            f"{purkinje}.rootStimulus.startTime": "0.0",
            f"{purkinje}.rootStimulus.duration": "0.0",
            f"{purkinje}.rootStimulus.intensity": "0.0",
            f"{purkinje}.outputVariables.export": "(activationTime)",
            f"{coupling}.conductionNetworkDomain": cls._NETWORK_NAME,
            f"{coupling}.couplingMode": "unidirectional",
            f"{coupling}.electroDomainCoupler": "eikonalMonodomainPvjCoupler",
            # eikonalMonodomainPvjCoupler.C:47-55 FatalErrors at construction
            # when rPvj is absent -- no graph-file escape, unlike
            # reactionDiffusionPvjCoupler. Without it this fixture built a
            # dictionary the solver would refuse to start on.
            f"{coupling}.rPvj": "1e5",
            f"{purkinje}.conductionSystemSolver": conduction_system_solver,
        }

    def _build(self, conduction_system_solver: str) -> str:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties

        return build_electro_properties(
            {
                "myocardiumSolver": "monodomainSolver",
                "ionicModel": "BuenoOrovio",
                "tissue": "epicardialCells",
            },
            overrides=self._overrides(conduction_system_solver),
        )

    def test_restitution_specific_keys_appear_in_applicable_entries(self) -> None:
        """useEdgeConductance/referenceConductance must be selectable once
        conductionSystemSolver=restitutionEikonalSolver1D is set -- not
        silently dropped as inapplicable."""
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            resolve_context,
            select_applicable_entries,
        )

        context = resolve_context(
            {"myocardiumSolver": "monodomainSolver", "ionicModel": "BuenoOrovio",
             "tissue": "epicardialCells"},
            overrides=self._overrides("restitutionEikonalSolver1D"),
        )
        paths = {e.driver_path for e in select_applicable_entries(context)}
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>."
            "purkinjeGraphModelCoeffs.useEdgeConductance",
            paths,
        )
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>."
            "purkinjeGraphModelCoeffs.referenceConductance",
            paths,
        )

    def test_restitution_build_emits_solver_specific_keys(self) -> None:
        """The synthesised dict must contain the restitution-only keys with
        their typical_value fallback, keyed under the concrete instance
        name (not the "<name>" template)."""
        text = self._build("restitutionEikonalSolver1D")
        self.assertIn("conductionSystemSolver restitutionEikonalSolver1D;", text)
        self.assertIn("useEdgeConductance true;", text)
        self.assertIn("referenceConductance 1.0;", text)

    def test_restitution_build_differs_from_plain_eikonal_by_more_than_the_solver_name(self) -> None:
        """The two builds must NOT be identical apart from the solver-name
        line -- that was exactly the reported defect (a dict claiming to be
        restitutionEikonalSolver1D while carrying none of what makes it
        one)."""
        restitution_lines = self._build("restitutionEikonalSolver1D").splitlines()
        plain_lines = self._build("eikonalSolver1D").splitlines()

        differing = sum(
            1 for a, b in zip(restitution_lines, plain_lines) if a != b
        )
        self.assertGreater(
            differing, 1,
            "restitutionEikonalSolver1D and eikonalSolver1D builds differ "
            f"by only {differing} line(s) -- restitution-specific keys are "
            "being dropped again",
        )
        # useEdgeConductance/referenceConductance must be present in the
        # restitution build and absent from the plain eikonalSolver1D build.
        restitution_text = "\n".join(restitution_lines)
        plain_text = "\n".join(plain_lines)
        self.assertIn("useEdgeConductance", restitution_text)
        self.assertNotIn("useEdgeConductance", plain_text)
        self.assertIn("referenceConductance", restitution_text)
        self.assertNotIn("referenceConductance", plain_text)

    def test_plain_eikonalSolver1D_does_not_gain_restitution_keys(self) -> None:
        """Non-regression: plain eikonalSolver1D must still NOT carry the
        restitution-only keys (they must stay conditional, not become
        unconditionally applicable)."""
        text = self._build("eikonalSolver1D")
        self.assertIn("conductionSystemSolver eikonalSolver1D;", text)
        self.assertNotIn("useEdgeConductance", text)
        self.assertNotIn("referenceConductance", text)


class TestRegenerateElectroProperties(unittest.TestCase):
    """`regenerate_electro_properties` -- the "route myocardiumSolver
    through the override channel" mechanism: rewrite an existing
    electroProperties file in place with one selector switched, rather
    than key-patching (which cannot rename the active <solver>Coeffs
    sub-block or drop now-illegal siblings)."""

    def test_module_exposes_regenerate_electro_properties(self) -> None:
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            regenerate_electro_properties,
        )
        self.assertTrue(callable(regenerate_electro_properties))

    def test_rejects_a_non_selector_driver_path(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            regenerate_electro_properties,
        )
        with tempfile.TemporaryDirectory() as d:
            p = self._write(d, build_electro_properties(
                {"myocardiumSolver": "singleCellSolver", "ionicModel": "AlievPanfilov", "tissue": "myocyte"},
            ))
            with self.assertRaises(ValueError) as exc:
                regenerate_electro_properties(p, "$ELECTRO_MODEL_COEFFS.solutionAlgorithm", "implicit")
            self.assertIn("not a regeneration selector key", str(exc.exception))

    @staticmethod
    def _write(tmp_dir, text):
        p = Path(tmp_dir) / "electroProperties"
        p.write_text(text)
        return p

    def test_renames_the_active_coeffs_block(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            regenerate_electro_properties,
        )
        with tempfile.TemporaryDirectory() as d:
            p = self._write(d, build_electro_properties(
                {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
            ))
            regenerate_electro_properties(
                p, "myocardiumSolver", "singleCellSolver",
                {"$ELECTRO_MODEL_COEFFS.ionicModel": "TNNP", "$ELECTRO_MODEL_COEFFS.tissue": "epicardialCells"},
            )
            text = p.read_text()
            self.assertIn("myocardiumSolver singleCellSolver;", text)
            self.assertIn("singleCellSolverCoeffs", text)
            self.assertNotIn("monodomainSolverCoeffs", text)

    def test_prunes_a_selector_that_becomes_forbidden_under_the_new_value(self) -> None:
        """ionicModel/tissue are forbidden_when myocardiumSolver=eikonalSolver
        (eikonalSolverCoeffs has no such keys at all). A regeneration that
        did not drop them would hand build_electro_properties a
        self-contradictory context and (correctly) fail -- this is the
        pruning that makes the solver-switch case succeed."""
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            regenerate_electro_properties,
        )
        with tempfile.TemporaryDirectory() as d:
            p = self._write(d, build_electro_properties(
                {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
            ))
            regenerate_electro_properties(
                p, "myocardiumSolver", "eikonalSolver",
                {
                    "$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach": "true",
                    "$ELECTRO_MODEL_COEFFS.stimulusLocationMin": "(1e6 1e6 1e6)",
                    "$ELECTRO_MODEL_COEFFS.stimulusLocationMax": "(1e6 1e6 1e6)",
                },
            )
            text = p.read_text()
            self.assertIn("myocardiumSolver eikonalSolver;", text)
            self.assertIn("eikonalSolverCoeffs", text)
            self.assertNotIn("ionicModel", text)
            self.assertNotIn("tissue", text)

    def test_raises_when_new_solver_requires_a_key_with_no_default_and_no_extra_override(self) -> None:
        """eikonalSolver's stimulusLocationMin/Max and
        eikonalAdvectionDiffusionApproach are required with no
        typical_value -- a solver switch that doesn't supply them must
        fail loudly (build_electro_properties's own validator), not emit
        an invalid dict."""
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            regenerate_electro_properties,
        )
        with tempfile.TemporaryDirectory() as d:
            p = self._write(d, build_electro_properties(
                {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
            ))
            with self.assertRaises(ValueError) as exc:
                regenerate_electro_properties(p, "myocardiumSolver", "eikonalSolver")
            self.assertIn("stimulusLocationMin", str(exc.exception))

    _PURKINJE_OVERRIDES = {
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.conductionSystemDomain": "purkinjeGraphModel",
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs.conductionSystemSolver": "monodomain1DSolver",
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs.graphFile": "purkinjeGraph",
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs.ionicModel": "BuenoOrovio",
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs.tissue": "epicardialCells",
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs.rootStimulus.node": "0",
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs.rootStimulus.startTime": "0.01",
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs.rootStimulus.duration": "0.002",
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs.rootStimulus.intensity": "500000.0",
        "$ELECTRO_MODEL_COEFFS.domainCouplings.couplingA.electroDomainCoupler": "reactionDiffusionPvjCoupler",
        "$ELECTRO_MODEL_COEFFS.domainCouplings.couplingA.conductionNetworkDomain": "purkinjeNetwork",
        # pvjCoupler.C:90 reads couplingMode with a hard dict.get<word>, so
        # every coupler needs it; the fixture predates that being enforced.
        "$ELECTRO_MODEL_COEFFS.domainCouplings.couplingA.couplingMode": "unidirectional",
    }

    def test_carries_forward_a_dynamic_container_verbatim(self) -> None:
        """conductionNetworkDomains is dynamic_path=True. Even though
        parse_electro_properties can now structurally round-trip most of it
        into overrides, regenerate_electro_properties still carries it
        forward verbatim (see _capture_dynamic_containers) rather than
        resynthesising it -- this is the hazard originally flagged: a naive
        rebuild must not silently drop a case's Purkinje network.

        Switches to bidomainSolver, which -- like the monodomainSolver this
        file started with -- physically supports a monodomain1DSolver
        Purkinje network via reactionDiffusionPvjCoupler (see
        solver_coupling.SOLVER_COMPATIBILITY_RULES), so the carried-forward
        network remains valid post-switch and the regeneration succeeds."""
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            regenerate_electro_properties,
        )
        with tempfile.TemporaryDirectory() as d:
            p = self._write(d, build_electro_properties(
                {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
                overrides=self._PURKINJE_OVERRIDES,
            ))
            self.assertIn("conductionNetworkDomains", p.read_text())

            regenerate_electro_properties(
                p, "myocardiumSolver", "bidomainSolver",
                {"$ELECTRO_MODEL_COEFFS.ionicModel": "TNNP", "$ELECTRO_MODEL_COEFFS.tissue": "epicardialCells"},
            )
            text = p.read_text()
            self.assertIn("bidomainSolverCoeffs", text)
            # The Purkinje block must still be present, verbatim, under the
            # renamed coeffs scope -- not dropped.
            self.assertIn("conductionNetworkDomains", text)
            self.assertIn("purkinjeNetwork", text)
            self.assertIn("conductionSystemSolver monodomain1DSolver;", text)
            self.assertIn("graphFile", text)
            self.assertIn("purkinjeGraph", text)

    def test_switching_to_a_physically_incompatible_solver_raises(self) -> None:
        """singleCellSolver has no PDE domain at all -- SOLVER_COMPATIBILITY_RULES
        marks ANY Purkinje pairing under it unconditionally invalid. Carrying
        the same monodomain1DSolver Purkinje network forward into a
        singleCellSolver switch must now be rejected instead of silently
        shipping a physically-meaningless file: the final safety net in
        regenerate_electro_properties validates the actually-shipped file,
        after the carried-forward block is reinserted, not just the
        intermediate rebuild the carried-forward block is invisible to."""
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            build_electro_properties,
            regenerate_electro_properties,
        )
        with tempfile.TemporaryDirectory() as d:
            p = self._write(d, build_electro_properties(
                {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
                overrides=self._PURKINJE_OVERRIDES,
            ))
            with self.assertRaises(ValueError) as exc:
                regenerate_electro_properties(
                    p, "myocardiumSolver", "singleCellSolver",
                    {"$ELECTRO_MODEL_COEFFS.ionicModel": "TNNP", "$ELECTRO_MODEL_COEFFS.tissue": "epicardialCells"},
                )
            message = str(exc.exception)
            self.assertIn("singleCellSolver", message)
            self.assertIn("incompatible", message.lower())
            # The rejected attempt must leave the ORIGINAL file completely
            # untouched -- the rewrite is staged on a scratch copy and only
            # committed to `path` once the final check passes, so a
            # rejection here must never ship a known-invalid intermediate
            # result.
            from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
            self.assertEqual(
                parse_electro_properties(p)["selectors"]["myocardiumSolver"],
                "monodomainSolver",
            )

    def test_purkinje_niederer_monodomain_to_eikonal_end_to_end(self) -> None:
        """Acceptance scenario, on a throwaway copy of the real
        purkinjeNiedererEtAl2011 tutorial's monodomain fixture: a bare
        myocardiumSolver switch to eikonalSolver, carrying the existing
        Purkinje network (conductionSystemSolver=monodomain1DSolver,
        electroDomainCoupler=reactionDiffusionPvjCoupler) forward verbatim,
        must now be REJECTED.

        This is not a hypothetical: the real, hand-authored, committed
        electroProperties.eikonal fixture for this same tutorial does NOT
        carry that Purkinje network forward -- it uses a different, physically
        compatible pairing instead (conductionSystemSolver=eikonalSolver1D,
        electroDomainCoupler=eikonalPvjCoupler; see
        solver_coupling.SOLVER_COMPATIBILITY_RULES, which marks
        eikonalSolver+monodomain1DSolver invalid: "eikonal myocardium cannot
        couple to reaction-diffusion Purkinje"). Producing a genuinely valid
        eikonal file from this tutorial requires resynthesising the Purkinje
        network's own selector/coupler -- a separate, deliberate step left to
        ordinary $ELECTRO_MODEL_COEFFS.* overrides applied after this one
        (regenerate_electro_properties's own documented contract), not
        something a bare myocardiumSolver switch can produce by carrying the
        old network forward unmodified. This test previously asserted the
        opposite (that carrying it forward verbatim succeeds and its static
        leaves are acceptable) without ever checking whether the carried-
        forward result was physically valid; it wasn't."""
        import shutil
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.dict_builder import (
            parse_electro_properties,
            regenerate_electro_properties,
        )

        if not PURKINJE_ELECTRO_PROPERTIES_MONODOMAIN.exists():
            self.skipTest("tutorial fixture not present in this checkout")

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "electroProperties"
            shutil.copyfile(PURKINJE_ELECTRO_PROPERTIES_MONODOMAIN, p)

            with self.assertRaises(ValueError) as exc:
                regenerate_electro_properties(
                    p, "myocardiumSolver", "eikonalSolver",
                    {
                        "$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach": "true",
                        "$ELECTRO_MODEL_COEFFS.stimulusLocationMin": "(1e6 1e6 1e6)",
                        "$ELECTRO_MODEL_COEFFS.stimulusLocationMax": "(1e6 1e6 1e6)",
                    },
                )
            message = str(exc.exception)
            self.assertIn("eikonal", message.lower())
            self.assertIn("reaction-diffusion", message.lower())

            # Original file untouched: still the monodomain tutorial fixture.
            self.assertEqual(
                parse_electro_properties(p)["selectors"]["myocardiumSolver"],
                "monodomainSolver",
            )


if __name__ == "__main__":
    unittest.main()
