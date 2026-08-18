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
#     test_postprocess_phase
#
# Description
#     Tests the post-DAG hand-off stub and its JSON serialization.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the post-DAG hand-off.

`run_postprocess_phase` (single-run) and `run_postprocessing_module`
(sweep) currently do no real analysis -- they prove the hand-off point
exists and returns a stable, serializable shape. `build_sweep_context` is
the "brain": it reads sweep_manifest.json and verifies it against what is
actually on disk per case, so `run_postprocessing_module` never has to
re-derive "what ran and where" itself.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.postprocess_phase import (
    CaseRecord,
    PostprocessOutcome,
    SweepContext,
    build_sweep_context,
    find_postprocess_script,
    read_case_output_file,
    read_case_workflow_state,
    run_postprocess_phase,
    run_postprocessing_module,
)
from openfoam_driver.core.runtime.sweep_manifest import (
    CaseManifestEntry,
    SweepManifest,
    write_manifest,
)


class PostprocessPhaseTests(unittest.TestCase):
    def test_returns_stub_outcome(self) -> None:
        outcome = run_postprocess_phase(entry="singleCell", output_dir=Path("/tmp/out"))
        self.assertIsInstance(outcome, PostprocessOutcome)
        self.assertEqual(outcome.status, "stub")
        self.assertIn("singleCell", outcome.message)
        self.assertIn("/tmp/out", outcome.message)

    def test_to_json_round_trips_status_and_message(self) -> None:
        outcome = run_postprocess_phase(entry="niederer2012", output_dir=Path("/tmp/x"))
        payload = outcome.to_json()
        self.assertEqual(payload, {"status": outcome.status, "message": outcome.message})


def _write_sweep_manifest(output_dir: Path, cases: list[CaseManifestEntry]) -> None:
    manifest = SweepManifest(
        schema_version="1.0",
        sweep_spec_hash="sha256:deadbeef",
        created_at="2026-08-18T00:00:00+00:00",
        updated_at="2026-08-18T00:05:00+00:00",
        cases=cases,
    )
    write_manifest(output_dir / "sweep_manifest.json", manifest)


class BuildSweepContextTests(unittest.TestCase):
    def test_resolves_relative_workflow_state_path_and_lists_real_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            case_dir = output_dir / "case_a" / "postProcessing"
            case_dir.mkdir(parents=True)
            (case_dir / "workflow_state.json").write_text("{}")
            (case_dir / "activationTime.csv").write_text("t,v\n0,0\n")

            _write_sweep_manifest(output_dir, [
                CaseManifestEntry(
                    case_id="case_a",
                    resolved_axis_values={"dx_mm": 0.5},
                    override_hash="sha256:x",
                    run_document_path="case_a/run_document.json",
                    workflow_state_path="case_a/postProcessing/workflow_state.json",
                    status="completed",
                    outcome="fresh",
                    started_at="2026-08-18T00:00:00+00:00",
                    updated_at="2026-08-18T00:01:00+00:00",
                ),
            ])

            context = build_sweep_context(output_dir)

            self.assertIsInstance(context, SweepContext)
            self.assertEqual(context.case_count, 1)
            self.assertEqual(context.completed_count, 1)
            self.assertEqual(context.failed_count, 0)
            self.assertEqual(len(context.cases), 1)
            case = context.cases[0]
            self.assertIsInstance(case, CaseRecord)
            self.assertEqual(case.case_id, "case_a")
            self.assertEqual(case.resolved_axis_values, {"dx_mm": 0.5})
            self.assertEqual(case.case_output_dir, str(case_dir))
            self.assertEqual(
                case.output_files,
                ("activationTime.csv", "workflow_state.json"),
            )

    def test_resolves_setup_root_from_real_run_document(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            case_dir = output_dir / "case_a" / "postProcessing"
            case_dir.mkdir(parents=True)
            (case_dir / "workflow_state.json").write_text("{}")

            run_document_dir = output_dir / "case_a"
            (run_document_dir / "run_document.json").write_text(json.dumps({
                "launch": {"setupRoot": "/tutorials/someTutorial/setup"},
            }))

            _write_sweep_manifest(output_dir, [
                CaseManifestEntry(
                    case_id="case_a",
                    resolved_axis_values={},
                    override_hash="sha256:x",
                    run_document_path="case_a/run_document.json",
                    workflow_state_path="case_a/postProcessing/workflow_state.json",
                    status="completed",
                    outcome="fresh",
                    started_at="2026-08-18T00:00:00+00:00",
                    updated_at="2026-08-18T00:01:00+00:00",
                ),
            ])

            context = build_sweep_context(output_dir)

            self.assertEqual(context.cases[0].setup_root, "/tutorials/someTutorial/setup")

    def test_missing_run_document_yields_no_setup_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            case_dir = output_dir / "case_a" / "postProcessing"
            case_dir.mkdir(parents=True)
            (case_dir / "workflow_state.json").write_text("{}")

            _write_sweep_manifest(output_dir, [
                CaseManifestEntry(
                    case_id="case_a",
                    resolved_axis_values={},
                    override_hash="sha256:x",
                    run_document_path="case_a/run_document.json",  # never written
                    workflow_state_path="case_a/postProcessing/workflow_state.json",
                    status="completed",
                    outcome="fresh",
                    started_at="2026-08-18T00:00:00+00:00",
                    updated_at="2026-08-18T00:01:00+00:00",
                ),
            ])

            context = build_sweep_context(output_dir)

            self.assertIsNone(context.cases[0].setup_root)

    def test_resolves_absolute_workflow_state_path_outside_output_dir(self) -> None:
        # Entry-mode sweeps record an absolute workflow_state_path pointing
        # into the target tutorial's own case_root, unrelated to the sweep's
        # --output-dir tree (see sweep_runner._relative_or_absolute).
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output_dir = root / "sweep_out"
            output_dir.mkdir()
            tutorial_case_dir = root / "tutorials" / "someTutorial" / "outputs"
            tutorial_case_dir.mkdir(parents=True)
            (tutorial_case_dir / "workflow_state.json").write_text("{}")

            _write_sweep_manifest(output_dir, [
                CaseManifestEntry(
                    case_id="implicit_TNNP_DX0.5",
                    resolved_axis_values={"dx_values": [0.5]},
                    override_hash="sha256:y",
                    run_document_path="implicit_TNNP_DX0.5/run_document.json",
                    workflow_state_path=str(tutorial_case_dir / "workflow_state.json"),
                    status="completed",
                    outcome="fresh",
                    started_at="2026-08-18T00:00:00+00:00",
                    updated_at="2026-08-18T00:01:00+00:00",
                ),
            ])

            context = build_sweep_context(output_dir)

            case = context.cases[0]
            self.assertEqual(case.case_output_dir, str(tutorial_case_dir))
            self.assertEqual(case.output_files, ("workflow_state.json",))

    def test_missing_workflow_state_file_yields_no_output_dir_or_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            _write_sweep_manifest(output_dir, [
                CaseManifestEntry(
                    case_id="never_ran",
                    resolved_axis_values={},
                    override_hash="sha256:z",
                    run_document_path="never_ran/run_document.json",
                    workflow_state_path="never_ran/postProcessing/workflow_state.json",
                    status="failed",
                    outcome="fresh",
                    started_at="2026-08-18T00:00:00+00:00",
                    updated_at="2026-08-18T00:01:00+00:00",
                ),
            ])

            context = build_sweep_context(output_dir)

            case = context.cases[0]
            self.assertIsNone(case.case_output_dir)
            self.assertEqual(case.output_files, ())
            self.assertEqual(context.failed_count, 1)


class FindPostprocessScriptTests(unittest.TestCase):
    def test_finds_script_defining_run_postprocessing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            setup_root = Path(tmp)
            (setup_root / "table_summary.py").write_text(
                "def run_postprocessing(*, output_dir, **kwargs):\n    return []\n"
            )
            self.assertEqual(
                find_postprocess_script(setup_root),
                setup_root / "table_summary.py",
            )

    def test_returns_none_when_no_script_matches(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            setup_root = Path(tmp)
            (setup_root / "helpers.py").write_text("def not_it():\n    pass\n")
            self.assertIsNone(find_postprocess_script(setup_root))

    def test_returns_none_for_nonexistent_directory(self) -> None:
        self.assertIsNone(find_postprocess_script(Path("/no/such/setup/dir")))

    def test_picks_first_match_in_sorted_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            setup_root = Path(tmp)
            (setup_root / "z_script.py").write_text("def run_postprocessing(): pass\n")
            (setup_root / "a_script.py").write_text("def run_postprocessing(): pass\n")
            self.assertEqual(
                find_postprocess_script(setup_root),
                setup_root / "a_script.py",
            )

    def test_never_imports_the_candidate_file(self) -> None:
        # A syntax error would only surface on import/exec -- a pure text
        # search must not choke on it, since setup/ scripts are untrusted
        # tutorial content the brain should never execute just to look.
        with tempfile.TemporaryDirectory() as tmp:
            setup_root = Path(tmp)
            (setup_root / "broken.py").write_text(
                "def run_postprocessing(:\n    this is not valid python\n"
            )
            self.assertEqual(
                find_postprocess_script(setup_root),
                setup_root / "broken.py",
            )


class RunPostprocessingModuleTests(unittest.TestCase):
    def test_message_is_grounded_in_the_supplied_context_and_task(self) -> None:
        context = SweepContext(
            output_dir="/tmp/sweep_out",
            sweep_spec_hash="sha256:abc123",
            started_at="2026-08-18T00:00:00+00:00",
            finished_at="2026-08-18T00:05:00+00:00",
            case_count=1,
            completed_count=1,
            failed_count=0,
            cases=(
                CaseRecord(
                    case_id="case_a",
                    resolved_axis_values={"dx_mm": 0.5},
                    status="completed",
                    outcome="fresh",
                    workflow_state_path="/tmp/sweep_out/case_a/postProcessing/workflow_state.json",
                    case_output_dir="/tmp/sweep_out/case_a/postProcessing",
                    output_files=("activationTime.csv",),
                    setup_root=None,
                ),
            ),
        )

        outcome = run_postprocessing_module(context, task="check convergence")

        self.assertIsInstance(outcome, PostprocessOutcome)
        self.assertEqual(outcome.status, "stub")
        self.assertIn("sha256:abc123", outcome.message)
        self.assertIn("case_a", outcome.message)
        self.assertIn("1/1 completed", outcome.message)
        self.assertIn("check convergence", outcome.message)
        self.assertIn("no setup/ script", outcome.message)

    def test_discovers_case_setup_root_postprocess_script(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            setup_root = Path(tmp) / "setup"
            setup_root.mkdir()
            (setup_root / "table_summary.py").write_text(
                "def run_postprocessing(*, output_dir, setup_root=None, **kwargs):\n    return []\n"
            )

            context = SweepContext(
                output_dir="/tmp/sweep_out",
                sweep_spec_hash="sha256:def456",
                started_at="2026-08-18T00:00:00+00:00",
                finished_at="2026-08-18T00:05:00+00:00",
                case_count=1,
                completed_count=1,
                failed_count=0,
                cases=(
                    CaseRecord(
                        case_id="case_a",
                        resolved_axis_values={},
                        status="completed",
                        outcome="fresh",
                        workflow_state_path="/tmp/sweep_out/case_a/postProcessing/workflow_state.json",
                        case_output_dir="/tmp/sweep_out/case_a/postProcessing",
                        output_files=(),
                        setup_root=str(setup_root),
                    ),
                ),
            )

            outcome = run_postprocessing_module(context, task="summarize")

            self.assertIn("found table_summary.py", outcome.message)


class OnDemandQueryTests(unittest.TestCase):
    def _build_real_context(self, output_dir: Path) -> SweepContext:
        case_dir = output_dir / "case_a" / "postProcessing"
        case_dir.mkdir(parents=True)
        (case_dir / "workflow_state.json").write_text(json.dumps({
            "status": "completed",
            "steps": [{"step_id": "solve", "status": "completed", "produced_artifacts": ["vm_series"]}],
        }))
        (case_dir / "activationTime.csv").write_text("t,v\n0,0\n")

        _write_sweep_manifest(output_dir, [
            CaseManifestEntry(
                case_id="case_a",
                resolved_axis_values={"dx_mm": 0.5},
                override_hash="sha256:x",
                run_document_path="case_a/run_document.json",
                workflow_state_path="case_a/postProcessing/workflow_state.json",
                status="completed",
                outcome="fresh",
                started_at="2026-08-18T00:00:00+00:00",
                updated_at="2026-08-18T00:01:00+00:00",
            ),
        ])
        return build_sweep_context(output_dir)

    def test_read_case_workflow_state_returns_full_step_detail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            context = self._build_real_context(Path(tmp))
            state = read_case_workflow_state(context, "case_a")
            self.assertEqual(state["status"], "completed")
            self.assertEqual(state["steps"][0]["produced_artifacts"], ["vm_series"])

    def test_read_case_workflow_state_unknown_case_id_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            context = self._build_real_context(Path(tmp))
            with self.assertRaises(KeyError):
                read_case_workflow_state(context, "does_not_exist")

    def test_read_case_output_file_returns_content(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            context = self._build_real_context(Path(tmp))
            content = read_case_output_file(context, "case_a", "activationTime.csv")
            self.assertEqual(content, "t,v\n0,0\n")

    def test_read_case_output_file_rejects_unverified_path(self) -> None:
        # Even though this exact file exists on disk (planted after the brain
        # already scanned), the brain never verified it -- must not be
        # readable via the query function. This is the guarantee the split
        # exists for: no path the brain hasn't confirmed is ever accessible.
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            context = self._build_real_context(output_dir)
            (output_dir / "case_a" / "postProcessing" / "planted_after_scan.txt").write_text("sneaky")
            with self.assertRaises(KeyError):
                read_case_output_file(context, "case_a", "planted_after_scan.txt")

    def test_read_case_output_file_unknown_case_id_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            context = self._build_real_context(Path(tmp))
            with self.assertRaises(KeyError):
                read_case_output_file(context, "does_not_exist", "activationTime.csv")


if __name__ == "__main__":
    unittest.main()
