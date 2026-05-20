"""Predictor contract tests (plan v2 phase 3a).

predict_data_artifacts is the single agent-facing answer to "what raw data
will/did this run produce?". It must:

* derive artifacts by composing existing catalogs (ionic_model_catalog,
  active_tension_catalog, dict_entries) — no solver-aware branching that
  re-implements catalog logic;
* merge a tutorial-supplied static override hook
  (``spec.metadata['expected_artifacts']``);
* never raise on unknown solvers or missing files — agents call it in
  exploratory contexts before a run has produced anything.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.artifacts import predict_data_artifacts
from openfoam_driver.core.runtime.models import (
    CaseConfig,
    DataArtifact,
    TutorialSpec,
    expand_path_pattern,
)


def _make_spec(
    case_root: Path,
    *,
    expected_artifacts: tuple[DataArtifact, ...] = (),
    cases: tuple[CaseConfig, ...] = (CaseConfig("only", {}),),
) -> TutorialSpec:
    metadata: dict[str, object] = {}
    if expected_artifacts:
        metadata["expected_artifacts"] = expected_artifacts
    return TutorialSpec(
        name="fixture",
        case_root=case_root,
        setup_root=case_root,
        output_dir=case_root / "output",
        build_cases=lambda: list(cases),
        apply_case=lambda _c, _case: None,
        run_case=lambda _c, _s, _case: None,
        metadata=metadata,
    )


def _write_single_cell_electro_properties(
    case_root: Path,
    *,
    ionic_model: str = "AlievPanfilov",
    tissue: str = "myocyte",
) -> None:
    """Synthesize the minimum-viable single-cell electroProperties.

    Mirrors the structure of
    tutorials/singleCellprotocols/singleCell/constant/electroProperties
    closely enough for the line-based parsers used by the rest of the
    driver to find the relevant keys.
    """
    constant = case_root / "constant"
    constant.mkdir(parents=True, exist_ok=True)
    (constant / "electroProperties").write_text(
        "myocardiumSolver singleCellSolver;\n"
        "singleCellSolverCoeffs\n"
        "{\n"
        f"    ionicModel    {ionic_model};\n"
        f"    tissue        {tissue};\n"
        "    solutionAlgorithm explicit;\n"
        "}\n"
    )


class TestPredictorSingleCell(unittest.TestCase):
    def test_emits_artifact_with_variables_from_catalog(self) -> None:
        """AlievPanfilov advertises states ('u', 'recovery_r') in the ionic
        model catalog — the predictor must source variables from there
        rather than redefining them locally."""
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_single_cell_electro_properties(
                case_root, ionic_model="AlievPanfilov"
            )
            spec = _make_spec(case_root)

            artifacts = predict_data_artifacts(case_root, spec)
            self.assertEqual(len(artifacts), 1)
            (artifact,) = artifacts
            self.assertEqual(artifact.produced_by, "singleCellSolver")
            self.assertIn("u", artifact.variables)
            self.assertIn("recovery_r", artifact.variables)

    def test_variables_change_with_ionic_model(self) -> None:
        """Convergence guard: the predictor must return different variables
        for different ionic models. Post-§3d-1 the variables come from
        recommended_exports (or the export list when declared); the exact
        names differ across models, proving catalog consultation."""
        with tempfile.TemporaryDirectory() as temp_a, tempfile.TemporaryDirectory() as temp_b:
            case_a = Path(temp_a) / "case"
            case_b = Path(temp_b) / "case"
            case_a.mkdir()
            case_b.mkdir()
            _write_single_cell_electro_properties(case_a, ionic_model="TNNP")
            _write_single_cell_electro_properties(case_b, ionic_model="AlievPanfilov")

            (artifact_tnnp,) = predict_data_artifacts(case_a, _make_spec(case_a))
            (artifact_ap,) = predict_data_artifacts(case_b, _make_spec(case_b))

            self.assertNotEqual(
                set(artifact_tnnp.variables), set(artifact_ap.variables),
                "predictor returned identical variables for two different "
                "ionic models — catalog consultation is broken",
            )
            # TNNP.recommended_exports references a calcium variable;
            # AlievPanfilov has no calcium.
            self.assertIn("calcium_Cai", artifact_tnnp.variables)

    def test_unknown_ionic_model_returns_empty_variables(self) -> None:
        """The predictor must not raise on a model name absent from the
        catalog — agents may mutate dicts to an as-yet-undefined model."""
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_single_cell_electro_properties(
                case_root, ionic_model="NotARealModel"
            )
            spec = _make_spec(case_root)

            artifacts = predict_data_artifacts(case_root, spec)
            self.assertEqual(len(artifacts), 1)
            (artifact,) = artifacts
            self.assertEqual(artifact.variables, ())


class TestPredictorMergesStaticOverride(unittest.TestCase):
    def test_static_expected_artifacts_passed_through(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_single_cell_electro_properties(case_root)
            static = DataArtifact(
                artifact_id="exact_error_norm",
                path_pattern="postProcessing/exact_error.json",
                format="json_summary",
                description="L2 error vs analytic solution",
            )
            spec = _make_spec(case_root, expected_artifacts=(static,))

            artifacts = predict_data_artifacts(case_root, spec)
            ids = {a.artifact_id for a in artifacts}
            self.assertIn("exact_error_norm", ids)

    def test_static_wins_on_artifact_id_collision(self) -> None:
        """spec.metadata['expected_artifacts'] is the authoring escape hatch.
        When a static entry shares an artifact_id with a derived one, the
        static description must win — the human knew something the predictor
        could not derive."""
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_single_cell_electro_properties(case_root)

            # Discover what the derived artifact_id is, then collide on it.
            derived = predict_data_artifacts(case_root, _make_spec(case_root))
            self.assertEqual(len(derived), 1)
            colliding_id = derived[0].artifact_id

            static = DataArtifact(
                artifact_id=colliding_id,
                path_pattern="custom/path",
                format="csv_probe",
                description="hand-authored override",
            )
            spec = _make_spec(case_root, expected_artifacts=(static,))

            artifacts = predict_data_artifacts(case_root, spec)
            by_id = {a.artifact_id: a for a in artifacts}
            self.assertEqual(by_id[colliding_id].description, "hand-authored override")
            self.assertEqual(by_id[colliding_id].path_pattern, "custom/path")


def _write_pde_electro_properties(
    case_root: Path,
    *,
    solver: str,
    ionic_model: str = "TNNP",
    export_list: tuple[str, ...] | None = None,
) -> None:
    """Synthesize a monodomain/bidomain electroProperties shell.

    Mirrors tutorials/manufacturedSolutions/{monodomain,bidomain}/constant/electroProperties
    closely enough for the line-based parsers; everything not relevant to the
    predictor is omitted. When ``export_list`` is supplied, an
    ``outputVariables.ionic.export ( ... )`` block is injected to exercise
    the §3d-1 filtering path.
    """
    (case_root / "constant").mkdir(parents=True, exist_ok=True)
    body = (
        f"myocardiumSolver  {solver};\n"
        f"{solver}Coeffs\n"
        "{\n"
        f"    ionicModel    {ionic_model};\n"
        "    solutionAlgorithm implicit;\n"
    )
    if export_list is not None:
        body += (
            "    outputVariables\n"
            "    {\n"
            "        ionic\n"
            "        {\n"
            f"            export ({' '.join(export_list)});\n"
            "        }\n"
            "    }\n"
        )
    body += "}\n"
    (case_root / "constant" / "electroProperties").write_text(body)


def _write_eikonal_electro_properties(case_root: Path) -> None:
    """Eikonal cases do not declare an ionicModel (constraint enforced by
    dict_entries.py); the predictor must not require one."""
    (case_root / "constant").mkdir(parents=True, exist_ok=True)
    (case_root / "constant" / "electroProperties").write_text(
        "myocardiumSolver eikonalSolver;\n"
        "eikonalSolverCoeffs\n"
        "{\n"
        "    conductivity (1 0 0  0 1 0  0 0 1);\n"
        "}\n"
    )


class TestPredictorMonodomain(unittest.TestCase):
    def test_emits_time_series_with_catalog_variables(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_pde_electro_properties(
                case_root, solver="monodomainSolver", ionic_model="TNNP"
            )
            spec = _make_spec(case_root)

            artifacts = predict_data_artifacts(case_root, spec)
            self.assertEqual(len(artifacts), 1)
            (artifact,) = artifacts
            self.assertEqual(artifact.produced_by, "monodomainSolver")
            self.assertEqual(artifact.format, "openfoam_time_dirs")
            self.assertTrue(artifact.time_indexed)
            self.assertIn("Vm", artifact.variables)
            # Catalog-sourced: TNNP.recommended_exports has calcium_Cai;
            # this proves the variables come from the catalog rather than
            # hard-coded in the handler.
            self.assertIn("calcium_Cai", artifact.variables)

    def test_pattern_uses_time_placeholder(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_pde_electro_properties(
                case_root, solver="monodomainSolver", ionic_model="AlievPanfilov"
            )
            spec = _make_spec(case_root)

            (artifact,) = predict_data_artifacts(case_root, spec)
            self.assertIn("{time}", artifact.path_pattern)


class TestPredictorBidomain(unittest.TestCase):
    def test_emits_phi_e_and_phi_i_in_variables(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_pde_electro_properties(
                case_root, solver="bidomainSolver", ionic_model="TNNP"
            )
            spec = _make_spec(case_root)

            (artifact,) = predict_data_artifacts(case_root, spec)
            self.assertEqual(artifact.produced_by, "bidomainSolver")
            self.assertIn("Vm", artifact.variables)
            self.assertIn("phiE", artifact.variables)
            self.assertIn("phiI", artifact.variables)


class TestPredictorEikonal(unittest.TestCase):
    def test_emits_psi_field_without_ionic_model_lookup(self) -> None:
        """Eikonal cases do not declare ionicModel — the predictor must
        produce a sensible artifact regardless (no KeyError fallthrough)."""
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_eikonal_electro_properties(case_root)
            spec = _make_spec(case_root)

            (artifact,) = predict_data_artifacts(case_root, spec)
            self.assertEqual(artifact.produced_by, "eikonalSolver")
            self.assertIn("psi", artifact.variables)
            self.assertIn("Vm", artifact.variables)


class TestPredictorExportListFiltering(unittest.TestCase):
    """Plan §3d-1: predictor must report what will actually be on disk, not
    the catalog superset. When ``outputVariables.ionic.export`` is declared,
    the exported subset wins. When absent, the catalog's
    ``recommended_exports`` is the fallback."""

    def test_export_list_overrides_catalog_states_for_monodomain(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_pde_electro_properties(
                case_root,
                solver="monodomainSolver",
                ionic_model="monodomainFDAManufactured",
                export_list=("u1", "u2", "u3"),
            )
            spec = _make_spec(case_root)
            (artifact,) = predict_data_artifacts(case_root, spec)
            # Solver-provided Vm is always present; ionic part filtered.
            self.assertEqual(artifact.variables, ("Vm", "u1", "u2", "u3"))

    def test_export_list_overrides_catalog_states_for_bidomain(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_pde_electro_properties(
                case_root,
                solver="bidomainSolver",
                ionic_model="TNNP",
                export_list=("V", "Cai"),
            )
            spec = _make_spec(case_root)
            (artifact,) = predict_data_artifacts(case_root, spec)
            self.assertEqual(
                artifact.variables, ("Vm", "phiE", "phiI", "V", "Cai"),
            )

    def test_export_list_overrides_catalog_for_single_cell(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            (case_root / "constant").mkdir()
            (case_root / "constant" / "electroProperties").write_text(
                "myocardiumSolver singleCellSolver;\n"
                "singleCellSolverCoeffs\n"
                "{\n"
                "    ionicModel    AlievPanfilov;\n"
                "    outputVariables\n"
                "    {\n"
                "        ionic\n"
                "        {\n"
                "            export (Vm s);\n"
                "        }\n"
                "    }\n"
                "}\n"
            )
            spec = _make_spec(case_root)
            (artifact,) = predict_data_artifacts(case_root, spec)
            # singleCell has no solver-provided PDE prefix; export list is
            # the entire variable set.
            self.assertEqual(artifact.variables, ("Vm", "s"))

    def test_missing_export_list_falls_back_to_recommended_exports(self) -> None:
        """AlievPanfilov.recommended_exports = ('u', 'recovery_r') in
        ionic_model_catalog.py. With no export declaration, the predictor
        must use the catalog fallback instead of states + algebraic."""
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_single_cell_electro_properties(
                case_root, ionic_model="AlievPanfilov"
            )
            spec = _make_spec(case_root)
            (artifact,) = predict_data_artifacts(case_root, spec)
            self.assertEqual(artifact.variables, ("u", "recovery_r"))


class TestPredictorManufacturedFdaRoundTrip(unittest.TestCase):
    """Plan §3.2 fixture: a verification tutorial declares analytic-error
    artifacts statically while the predictor still derives the field series.
    Demonstrates that solver-derived + spec-declared artifacts coexist."""

    def test_both_derived_and_static_present(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            _write_pde_electro_properties(
                case_root,
                solver="monodomainSolver",
                ionic_model="monodomainFDAManufactured",
            )
            error_norm = DataArtifact(
                artifact_id="exact_error_norm",
                path_pattern="postProcessing/exact_error.json",
                format="json_summary",
                description="L2 error vs analytic manufactured solution",
                produced_by="manufacturedFDAMonodomainVerifier",
            )
            spec = _make_spec(case_root, expected_artifacts=(error_norm,))

            artifacts = predict_data_artifacts(case_root, spec)
            by_id = {a.artifact_id: a for a in artifacts}
            self.assertIn("myocardium_time_series", by_id)
            self.assertIn("exact_error_norm", by_id)
            # Derived artifact gets catalog variables for the manufactured model.
            derived = by_id["myocardium_time_series"]
            self.assertEqual(derived.produced_by, "monodomainSolver")
            self.assertIn("Vm", derived.variables)


class TestPredictorPathPatternContract(unittest.TestCase):
    """Every solver handler must emit patterns that round-trip through
    expand_path_pattern. This locks the predictor against drift from the
    closed placeholder enum in models.py."""

    def _all_solver_fixtures(self) -> list[tuple[str, callable]]:
        return [
            ("singleCellSolver",
             lambda case_root: _write_single_cell_electro_properties(case_root)),
            ("monodomainSolver",
             lambda case_root: _write_pde_electro_properties(
                 case_root, solver="monodomainSolver", ionic_model="TNNP")),
            ("bidomainSolver",
             lambda case_root: _write_pde_electro_properties(
                 case_root, solver="bidomainSolver", ionic_model="TNNP")),
            ("eikonalSolver",
             lambda case_root: _write_eikonal_electro_properties(case_root)),
        ]

    def test_every_handler_emits_expandable_patterns(self) -> None:
        for solver_name, writer in self._all_solver_fixtures():
            with self.subTest(solver=solver_name), tempfile.TemporaryDirectory() as temp:
                case_root = Path(temp) / "case"
                case_root.mkdir()
                writer(case_root)
                spec = _make_spec(case_root)
                artifacts = predict_data_artifacts(case_root, spec)
                self.assertGreater(
                    len(artifacts), 0,
                    f"{solver_name} produced no artifacts — fixture/handler mismatch",
                )
                for artifact in artifacts:
                    expanded = expand_path_pattern(
                        artifact.path_pattern,
                        case_id="probeCase",
                        time="0.001",
                    )
                    self.assertIsInstance(expanded, str)
                    self.assertNotIn(
                        "{", expanded,
                        f"{solver_name} pattern {artifact.path_pattern!r} expanded "
                        f"to {expanded!r} but still contains a brace",
                    )


class TestPredictorGracefulFallback(unittest.TestCase):
    def test_missing_electro_properties_returns_only_static(self) -> None:
        """Before a run mutates the case, electroProperties may not yet
        exist (or live under setup/). The predictor must degrade to the
        static hook rather than raising."""
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            static = DataArtifact(
                artifact_id="placeholder",
                path_pattern="postProcessing/placeholder.csv",
                format="csv_probe",
            )
            spec = _make_spec(case_root, expected_artifacts=(static,))

            artifacts = predict_data_artifacts(case_root, spec)
            self.assertEqual(len(artifacts), 1)
            self.assertEqual(artifacts[0].artifact_id, "placeholder")

    def test_unknown_solver_returns_only_static(self) -> None:
        """Unhandled solvers (e.g., a future addition) must not error.
        Phases will fill in coverage incrementally."""
        with tempfile.TemporaryDirectory() as temp:
            case_root = Path(temp) / "case"
            case_root.mkdir()
            (case_root / "constant").mkdir()
            (case_root / "constant" / "electroProperties").write_text(
                "myocardiumSolver totallyMadeUpSolver;\n"
                "totallyMadeUpSolverCoeffs { }\n"
            )
            spec = _make_spec(case_root)

            artifacts = predict_data_artifacts(case_root, spec)
            self.assertEqual(artifacts, ())


if __name__ == "__main__":
    unittest.main()
