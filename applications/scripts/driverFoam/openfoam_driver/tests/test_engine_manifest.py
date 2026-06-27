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
#     test_engine_manifest
#
# Description
#     Tests engine manifest logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import os
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

from openfoam_driver.core.runtime.engine import DriverEngine
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text())


class TestDriverEngineManifest(unittest.TestCase):
    def test_run_simulations_writes_rich_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            applied_cases: list[str] = []
            ran_cases: list[str] = []

            def build_cases() -> list[CaseConfig]:
                return [
                    CaseConfig("caseA", {"value": 1}),
                    CaseConfig("caseB", {"value": 2}),
                ]

            def apply_case(_: Path, case: CaseConfig) -> None:
                applied_cases.append(case.case_id)

            def run_case(_: Path, __: Path, case: CaseConfig) -> None:
                ran_cases.append(case.case_id)

            spec = TutorialSpec(
                name="dummy",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=build_cases,
                apply_case=apply_case,
                run_case=run_case,
            )

            engine = DriverEngine(spec=spec, requested_action="all")
            results = engine.run_simulations()

            manifest = _load_json(output_dir / "run_manifest.json")
            report_path = output_dir / "run_report.md"
            self.assertEqual(applied_cases, ["caseA", "caseB"])
            self.assertEqual(ran_cases, ["caseA", "caseB"])
            self.assertEqual(len(results), 2)
            self.assertEqual(manifest["schema_version"], "2.3")
            self.assertEqual(manifest["requested_action"], "all")
            self.assertEqual(manifest["entry"], "dummy")
            self.assertIsNone(manifest["entry_kind"])
            self.assertEqual(manifest["status"], "completed")
            self.assertEqual(manifest["postprocess_status"], "not_started")
            self.assertEqual(manifest["completed_cases"], 2)
            self.assertEqual(manifest["failed_cases"], 0)
            self.assertIsNotNone(manifest["run_id"])
            self.assertIsNotNone(manifest["started_at_utc"])
            self.assertIsNotNone(manifest["finished_at_utc"])
            self.assertEqual(manifest["human_report_path"], str(report_path))
            self.assertEqual(manifest["results"][0]["index"], 1)
            self.assertEqual(manifest["results"][1]["total_cases"], 2)
            self.assertTrue(report_path.exists())
            report_text = report_path.read_text()
            self.assertIn("# driverFoam Run Report", report_text)
            self.assertIn("- Entry: `dummy`", report_text)
            self.assertIn("### caseA", report_text)
            self.assertIn("### caseB", report_text)

    def test_run_postprocess_writes_manifest_and_plots_path(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            def postprocess(_: Path, out_dir: Path) -> None:
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / "plots.json").write_text('{"plots": []}')

            spec = TutorialSpec(
                name="dummy",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=lambda: [],
                apply_case=lambda _case_root, _case: None,
                run_case=lambda _case_root, _setup_root, _case: None,
                postprocess=postprocess,
            )

            engine = DriverEngine(spec=spec, requested_action="post")
            engine.run_postprocess()

            manifest = _load_json(output_dir / "run_manifest.json")
            self.assertEqual(manifest["requested_action"], "post")
            self.assertEqual(manifest["entry"], "dummy")
            self.assertEqual(manifest["status"], "completed")
            self.assertEqual(manifest["postprocess_status"], "completed")
            self.assertEqual(manifest["total_cases"], 0)
            self.assertEqual(manifest["results"], [])
            self.assertEqual(
                manifest["plots_manifest_path"],
                str(output_dir / "plots.json"),
            )
            self.assertTrue((output_dir / "run_report.md").exists())

    def test_run_postprocess_collects_outputs_before_postprocessing(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            collected: list[str] = []
            seen_in_postprocess: list[str] = []

            def collect_outputs(case_dir: Path, out_dir: Path) -> None:
                self.assertEqual(case_dir, case_root)
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / "3D_80_cells_implicit.dat").write_text("fresh-output")
                collected.append("ok")

            def postprocess(_: Path, out_dir: Path) -> None:
                seen_in_postprocess.extend(path.name for path in out_dir.glob("*.dat"))
                (out_dir / "plots.json").write_text('{"plots": []}')

            spec = TutorialSpec(
                name="dummy",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=lambda: [],
                apply_case=lambda _case_root, _case: None,
                run_case=lambda _case_root, _setup_root, _case: None,
                collect_outputs=collect_outputs,
                postprocess=postprocess,
            )

            engine = DriverEngine(spec=spec, requested_action="post")
            engine.run_postprocess()

            self.assertEqual(collected, ["ok"])
            self.assertEqual(seen_in_postprocess, ["3D_80_cells_implicit.dat"])
            self.assertTrue((output_dir / "3D_80_cells_implicit.dat").exists())

class TestDriverEngineManifestSchemaVersion(unittest.TestCase):
    """run_manifest.json schema version contract.

    v2.x is additive-only: any 2.x manifest must contain every key that a
    2.1 consumer expected to read. Bumps are reserved for new optional
    fields, never for renames or deletions.
    """

    # Locked v2.1 key set. Removing any key here is a breaking change and
    # requires a major version bump (3.0) plus a migration helper.
    _V2_1_REQUIRED_KEYS: frozenset[str] = frozenset({
        "schema_version", "run_id", "requested_action", "entry", "entry_kind",
        "entry_path", "source_type", "workflow_family", "case_root",
        "setup_root", "output_dir", "dry_run", "continue_on_error", "status",
        "postprocess_status", "current_case_id", "started_at_utc",
        "updated_at_utc", "finished_at_utc", "total_cases", "planned_cases",
        "completed_cases", "failed_cases", "error", "plots_manifest_path",
        "results", "human_report_path",
    })

    def test_manifest_schema_version_is_2_3_and_additive(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            spec = TutorialSpec(
                name="schema",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=lambda: [CaseConfig("only", {})],
                apply_case=lambda _c, _case: None,
                run_case=lambda _c, _s, _case: None,
            )
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            manifest = _load_json(output_dir / "run_manifest.json")
            self.assertEqual(manifest["schema_version"], "2.3")
            missing = self._V2_1_REQUIRED_KEYS - set(manifest)
            self.assertEqual(
                missing, set(),
                f"v2.2 manifest is missing v2.1 keys: {sorted(missing)}",
            )


class TestActionEventsJsonl(unittest.TestCase):
    """action_events.jsonl contract: one well-formed JSON object per line.

    Agents tail this file to track run progress in real time. The invariant
    is that each line is a single ``json.dumps(event) + "\n"`` write, so a
    line is either fully present or not present at all — never partial.
    """

    def test_every_line_is_standalone_json(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            spec = TutorialSpec(
                name="events",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=lambda: [
                    CaseConfig("e1", {}), CaseConfig("e2", {}), CaseConfig("e3", {}),
                ],
                apply_case=lambda _c, _case: None,
                run_case=lambda _c, _s, _case: None,
            )
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            events_path = output_dir / "action_events.jsonl"
            self.assertTrue(events_path.exists())
            lines = events_path.read_text().splitlines()
            self.assertGreater(len(lines), 0)

            events = []
            for idx, line in enumerate(lines):
                self.assertNotIn(
                    "\n", line,
                    f"line {idx} contains an embedded newline — events must be single-line",
                )
                try:
                    events.append(json.loads(line))
                except json.JSONDecodeError as exc:
                    self.fail(f"line {idx} is not valid JSON: {exc!r}\nline: {line!r}")

            # Every event carries the agent-required contract fields.
            for evt in events:
                self.assertIn("event", evt)
                self.assertIn("timestamp_utc", evt)
                self.assertIn("run_id", evt)

            event_types = [evt["event"] for evt in events]
            self.assertEqual(event_types[0], "sim_started")
            self.assertEqual(event_types[-1], "sim_finished")


class TestArtifactsManifest(unittest.TestCase):
    """artifacts_manifest.json contract (plan v2 §3b).

    Sidecar manifest written next to run_manifest.json, listing the predicted
    DataArtifacts for the current on-disk case state. Always written so
    agents have a stable polling target; the artifacts list may be empty.
    """

    def _write_single_cell_electro_properties(self, case_root: Path) -> None:
        (case_root / "constant").mkdir(parents=True, exist_ok=True)
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n"
            "{\n"
            "    ionicModel    AlievPanfilov;\n"
            "    tissue        myocyte;\n"
            "}\n"
        )

    def _build_spec(self, root: Path, *, with_electro: bool = True) -> TutorialSpec:
        case_root = root / "case"
        setup_root = root / "setup"
        output_dir = root / "output"
        case_root.mkdir()
        setup_root.mkdir()
        if with_electro:
            self._write_single_cell_electro_properties(case_root)
        return TutorialSpec(
            name="artifacts",
            case_root=case_root,
            setup_root=setup_root,
            output_dir=output_dir,
            build_cases=lambda: [CaseConfig("only", {})],
            apply_case=lambda _c, _case: None,
            run_case=lambda _c, _s, _case: None,
        )

    def test_artifacts_manifest_written_after_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            artifacts_path = spec.output_dir / "artifacts_manifest.json"
            self.assertTrue(
                artifacts_path.exists(),
                "engine must write artifacts_manifest.json as a sidecar to run_manifest.json",
            )
            payload = json.loads(artifacts_path.read_text())
            self.assertIn("schema_version", payload)
            self.assertIn("artifacts", payload)
            self.assertEqual(payload["run_id"],
                             json.loads((spec.output_dir / "run_manifest.json").read_text())["run_id"])

    def test_artifacts_match_predictor_output(self) -> None:
        from openfoam_driver.core.runtime.artifacts import predict_data_artifacts

        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            payload = json.loads(
                (spec.output_dir / "artifacts_manifest.json").read_text()
            )
            predicted = predict_data_artifacts(spec.case_root, spec)
            self.assertEqual(len(payload["artifacts"]), len(predicted))
            ids = {a["artifact_id"] for a in payload["artifacts"]}
            self.assertEqual(ids, {a.artifact_id for a in predicted})
            for emitted, expected in zip(payload["artifacts"], predicted):
                self.assertEqual(emitted["artifact_id"], expected.artifact_id)
                self.assertEqual(emitted["format"], expected.format)
                self.assertEqual(tuple(emitted["variables"]), expected.variables)
                self.assertEqual(emitted["produced_by"], expected.produced_by)

    def test_run_manifest_records_artifacts_path(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            run_manifest = json.loads(
                (spec.output_dir / "run_manifest.json").read_text()
            )
            self.assertEqual(
                run_manifest["artifacts_manifest_path"],
                str(spec.output_dir / "artifacts_manifest.json"),
            )

    def test_artifacts_manifest_written_atomically(self) -> None:
        """Mirrors the run_manifest atomicity contract — agents poll both."""
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            artifacts_path = spec.output_dir / "artifacts_manifest.json"

            tmp_existed_at_replace: list[bool] = []
            saw_artifacts_dst = False
            real_replace = os.replace

            def tracking_replace(src, dst):
                nonlocal saw_artifacts_dst
                src_path = Path(src)
                dst_path = Path(dst)
                if dst_path == artifacts_path:
                    saw_artifacts_dst = True
                    tmp_existed_at_replace.append(src_path.exists())
                    self.assertEqual(src_path.parent, artifacts_path.parent)
                    self.assertTrue(
                        src_path.name.startswith("artifacts_manifest.json")
                        and src_path.name.endswith(".tmp")
                    )
                return real_replace(src, dst)

            with mock.patch(
                "openfoam_driver.core.runtime.engine.os.replace",
                side_effect=tracking_replace,
            ):
                DriverEngine(spec=spec, requested_action="sim").run_simulations()

            self.assertTrue(saw_artifacts_dst,
                            "no os.replace targeting artifacts_manifest.json")
            self.assertTrue(all(tmp_existed_at_replace))
            leftovers = list(spec.output_dir.glob("artifacts_manifest.json*.tmp"))
            self.assertEqual(leftovers, [])

    def test_empty_artifacts_when_no_electro_properties(self) -> None:
        """File still written (stable polling target); artifacts list empty."""
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp), with_electro=False)
            DriverEngine(spec=spec, requested_action="sim").run_simulations()
            payload = json.loads(
                (spec.output_dir / "artifacts_manifest.json").read_text()
            )
            self.assertEqual(payload["artifacts"], [])


class TestDriverEngineManifestAtomicity(unittest.TestCase):
    """run_manifest.json must be readable by polling agents at any instant.

    The engine rewrites run_manifest.json many times during a sweep. A naive
    write_text() leaves a window where the file exists but contains a truncated
    or partial JSON document; concurrent readers raise json.JSONDecodeError.
    The contract: every write goes through a sibling .tmp file, then os.replace
    atomically swaps it into place.
    """

    def _build_spec(self, root: Path) -> TutorialSpec:
        case_root = root / "case"
        setup_root = root / "setup"
        output_dir = root / "output"
        case_root.mkdir()
        setup_root.mkdir()
        return TutorialSpec(
            name="atomic",
            case_root=case_root,
            setup_root=setup_root,
            output_dir=output_dir,
            build_cases=lambda: [
                CaseConfig("c1", {"v": 1}),
                CaseConfig("c2", {"v": 2}),
                CaseConfig("c3", {"v": 3}),
            ],
            apply_case=lambda _case_root, _case: None,
            run_case=lambda _case_root, _setup_root, _case: None,
        )

    def test_run_manifest_written_atomically_via_replace(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            spec = self._build_spec(Path(temp_dir))
            manifest_path = spec.output_dir / "run_manifest.json"

            replace_dsts: list[Path] = []
            tmp_existed_at_replace: list[bool] = []
            real_replace = os.replace

            def tracking_replace(src, dst):
                src_path = Path(src)
                dst_path = Path(dst)
                if dst_path == manifest_path:
                    replace_dsts.append(dst_path)
                    tmp_existed_at_replace.append(src_path.exists())
                    self.assertEqual(
                        src_path.parent,
                        manifest_path.parent,
                        "tmp file must be a sibling of the final manifest",
                    )
                    self.assertTrue(
                        src_path.name.startswith("run_manifest.json"),
                        f"tmp name should be derived from final name, got {src_path.name!r}",
                    )
                return real_replace(src, dst)

            engine = DriverEngine(spec=spec, requested_action="sim")
            with mock.patch(
                "openfoam_driver.core.runtime.engine.os.replace",
                side_effect=tracking_replace,
            ):
                engine.run_simulations()

            self.assertGreater(
                len(replace_dsts), 0,
                "expected at least one os.replace targeting run_manifest.json",
            )
            self.assertTrue(
                all(tmp_existed_at_replace),
                "tmp source must exist at the moment os.replace is called",
            )
            self.assertTrue(manifest_path.exists())
            # No leftover .tmp files after a successful run.
            leftover_tmps = list(spec.output_dir.glob("run_manifest.json*.tmp"))
            self.assertEqual(leftover_tmps, [], f"leftover tmp files: {leftover_tmps}")
            # Final document parses cleanly.
            _load_json(manifest_path)

    def test_concurrent_reader_never_sees_partial_manifest(self) -> None:
        """Agent-style poller reads run_manifest.json while engine writes it.

        Without atomicity, the reader will occasionally observe an empty or
        truncated file and raise JSONDecodeError. With atomic rename, every
        successful read returns a fully-formed manifest.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            setup_root = root / "setup"
            output_dir = root / "output"
            case_root.mkdir()
            setup_root.mkdir()

            # Slow the engine just enough that the poller observes
            # many manifest rewrites during the sweep.
            def slow_apply(_case_root: Path, _case: CaseConfig) -> None:
                time.sleep(0.002)

            spec = TutorialSpec(
                name="poller",
                case_root=case_root,
                setup_root=setup_root,
                output_dir=output_dir,
                build_cases=lambda: [
                    CaseConfig(f"c{i}", {"v": i}) for i in range(20)
                ],
                apply_case=slow_apply,
                run_case=lambda _case_root, _setup_root, _case: None,
            )

            manifest_path = output_dir / "run_manifest.json"
            stop_event = threading.Event()
            decode_errors: list[Exception] = []
            successful_reads = 0

            def poll_reader():
                nonlocal successful_reads
                while not stop_event.is_set():
                    if not manifest_path.exists():
                        continue
                    try:
                        json.loads(manifest_path.read_text())
                        successful_reads += 1
                    except json.JSONDecodeError as exc:
                        decode_errors.append(exc)
                    except FileNotFoundError:
                        # Tolerated only if manifest existed momentarily before
                        # the engine first wrote it; once written, atomic
                        # replace guarantees it never disappears.
                        pass

            reader = threading.Thread(target=poll_reader, daemon=True)
            reader.start()
            try:
                DriverEngine(spec=spec, requested_action="sim").run_simulations()
            finally:
                stop_event.set()
                reader.join(timeout=2.0)

            self.assertEqual(
                decode_errors, [],
                f"poller observed {len(decode_errors)} torn read(s); "
                f"first: {decode_errors[0] if decode_errors else None}",
            )
            self.assertGreater(
                successful_reads, 0,
                "poller never read the manifest — test scaffolding broken",
            )


class TestWorkflowDag(unittest.TestCase):
    """Tests for the workflow_dag field in run_manifest.json."""

    _EXAMPLE_DAG = {
        "steps": [
            {"id": "mesh", "command": "blockMesh", "depends_on": []},
            {"id": "solve", "command": "cardiacFoam", "depends_on": ["mesh"]},
        ]
    }

    def _build_spec_with_dag(self, root: Path) -> TutorialSpec:
        case_root = root / "case"
        setup_root = root / "setup"
        output_dir = root / "output"
        case_root.mkdir()
        setup_root.mkdir()
        return TutorialSpec(
            name="dag_test",
            case_root=case_root,
            setup_root=setup_root,
            output_dir=output_dir,
            build_cases=lambda: [CaseConfig("only", {})],
            apply_case=lambda _c, _case: None,
            run_case=lambda _c, _s, _case: None,
            metadata={"workflow_dag": self._EXAMPLE_DAG},
        )

    def _build_spec_without_dag(self, root: Path) -> TutorialSpec:
        case_root = root / "case"
        setup_root = root / "setup"
        output_dir = root / "output"
        case_root.mkdir()
        setup_root.mkdir()
        return TutorialSpec(
            name="no_dag",
            case_root=case_root,
            setup_root=setup_root,
            output_dir=output_dir,
            build_cases=lambda: [CaseConfig("only", {})],
            apply_case=lambda _c, _case: None,
            run_case=lambda _c, _s, _case: None,
        )

    def test_run_manifest_records_workflow_dag(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec_with_dag(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            manifest = _load_json(spec.output_dir / "run_manifest.json")
            self.assertIn("workflow_dag", manifest)
            self.assertEqual(manifest["workflow_dag"], self._EXAMPLE_DAG)
            # round-trip: steps are preserved exactly
            steps = manifest["workflow_dag"]["steps"]
            self.assertEqual(len(steps), 2)
            self.assertEqual(steps[0]["id"], "mesh")
            self.assertEqual(steps[0]["command"], "blockMesh")
            self.assertEqual(steps[0]["depends_on"], [])
            self.assertEqual(steps[1]["id"], "solve")
            self.assertEqual(steps[1]["command"], "cardiacFoam")
            self.assertEqual(steps[1]["depends_on"], ["mesh"])

    def test_run_manifest_workflow_dag_is_none_when_absent(self) -> None:
        """Spec without workflow_dag → manifest field is None (additive-only)."""
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec_without_dag(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            manifest = _load_json(spec.output_dir / "run_manifest.json")
            self.assertIn("workflow_dag", manifest)
            self.assertIsNone(manifest["workflow_dag"])

    def test_workflow_dag_does_not_break_v2_1_consumer(self) -> None:
        """A v2.2 manifest with workflow_dag still contains every v2.1 required key."""
        v2_1_required = TestDriverEngineManifestSchemaVersion._V2_1_REQUIRED_KEYS
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec_with_dag(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            manifest = _load_json(spec.output_dir / "run_manifest.json")
            self.assertEqual(manifest["schema_version"], "2.3")
            missing = v2_1_required - set(manifest)
            self.assertEqual(
                missing, set(),
                f"v2.2 manifest with workflow_dag is missing v2.1 keys: {sorted(missing)}",
            )


class TestArtifactsRealizedManifest(unittest.TestCase):
    """artifacts_realized.json contract (plan v2 §10).

    Sidecar written at terminal status reporting which predicted artifacts
    actually appeared on disk. Always present at the end of a non-dry run,
    so agents have a stable reconciliation target. Linked from
    run_manifest['artifacts_realized_path'].
    """

    def _write_single_cell_electro_properties(self, case_root: Path) -> None:
        (case_root / "constant").mkdir(parents=True, exist_ok=True)
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n"
            "{\n"
            "    ionicModel    AlievPanfilov;\n"
            "    tissue        myocyte;\n"
            "}\n"
        )

    def _build_spec(self, root: Path) -> TutorialSpec:
        case_root = root / "case"
        setup_root = root / "setup"
        output_dir = root / "output"
        case_root.mkdir()
        setup_root.mkdir()
        self._write_single_cell_electro_properties(case_root)
        return TutorialSpec(
            name="realized",
            case_root=case_root,
            setup_root=setup_root,
            output_dir=output_dir,
            build_cases=lambda: [CaseConfig("only", {})],
            apply_case=lambda _c, _case: None,
            run_case=lambda _c, _s, _case: None,
        )

    def _build_sweep_spec(self, root: Path) -> TutorialSpec:
        """Two-case sweep variant of _build_spec."""
        case_root = root / "case"
        setup_root = root / "setup"
        output_dir = root / "output"
        case_root.mkdir()
        setup_root.mkdir()
        self._write_single_cell_electro_properties(case_root)
        return TutorialSpec(
            name="sweep",
            case_root=case_root,
            setup_root=setup_root,
            output_dir=output_dir,
            build_cases=lambda: [CaseConfig("c1", {}), CaseConfig("c2", {})],
            apply_case=lambda _c, _case: None,
            run_case=lambda _c, _s, _case: None,
        )

    def test_artifacts_realized_written_at_terminal_status(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            realized = spec.output_dir / "artifacts_realized.json"
            self.assertTrue(realized.exists(), "artifacts_realized.json missing")
            payload = json.loads(realized.read_text())
            self.assertEqual(payload["schema_version"], "1.1")
            self.assertIn("cases", payload)
            self.assertGreaterEqual(len(payload["cases"]), 1)
            first = payload["cases"][0]
            self.assertIn("predicted_count", first)
            self.assertIn("matched_count", first)
            self.assertIn("missing_count", first)
            self.assertIn("artifacts", first)

    def test_realized_manifest_has_one_entry_per_case_in_sweep(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_sweep_spec(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            realized = json.loads(
                (spec.output_dir / "artifacts_realized.json").read_text()
            )
            self.assertEqual(realized["schema_version"], "1.1")
            self.assertIn("cases", realized)
            case_ids = [c["case_id"] for c in realized["cases"]]
            self.assertEqual(case_ids, ["c1", "c2"])

    def test_single_case_run_still_emits_cases_array(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            realized = json.loads(
                (spec.output_dir / "artifacts_realized.json").read_text()
            )
            self.assertEqual(realized["schema_version"], "1.1")
            self.assertIn("cases", realized)
            self.assertEqual(len(realized["cases"]), 1)
            self.assertEqual(realized["cases"][0]["case_id"], "only")

    def test_run_manifest_records_artifacts_realized_path(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            DriverEngine(spec=spec, requested_action="sim").run_simulations()

            run_manifest = json.loads(
                (spec.output_dir / "run_manifest.json").read_text()
            )
            self.assertEqual(
                run_manifest["artifacts_realized_path"],
                str(spec.output_dir / "artifacts_realized.json"),
            )

    def test_realized_manifest_is_written_atomically(self) -> None:
        """Same atomicity contract as run_manifest + artifacts_manifest —
        the realized manifest goes through .tmp + os.replace."""
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            realized_path = spec.output_dir / "artifacts_realized.json"

            tmp_existed_at_replace: list[bool] = []
            saw_realized_dst = False
            real_replace = os.replace

            def tracking_replace(src, dst):
                nonlocal saw_realized_dst
                src_path = Path(src)
                dst_path = Path(dst)
                if dst_path == realized_path:
                    saw_realized_dst = True
                    tmp_existed_at_replace.append(src_path.exists())
                    self.assertEqual(src_path.parent, realized_path.parent)
                    self.assertTrue(
                        src_path.name.startswith("artifacts_realized.json")
                        and src_path.name.endswith(".tmp"),
                        f"unexpected tmp name: {src_path.name}",
                    )
                return real_replace(src, dst)

            with mock.patch(
                "openfoam_driver.core.runtime.engine.os.replace",
                side_effect=tracking_replace,
            ):
                DriverEngine(spec=spec, requested_action="sim").run_simulations()

            self.assertTrue(saw_realized_dst)
            self.assertTrue(all(tmp_existed_at_replace))
            leftovers = list(spec.output_dir.glob("artifacts_realized.json*.tmp"))
            self.assertEqual(leftovers, [])

    def test_realized_manifest_skipped_on_dry_run(self) -> None:
        """Dry runs never execute case work — no realized artifacts to
        report. The sidecar must NOT be written."""
        with tempfile.TemporaryDirectory() as temp:
            spec = self._build_spec(Path(temp))
            DriverEngine(
                spec=spec, dry_run=True, requested_action="sim",
            ).run_simulations()
            realized = spec.output_dir / "artifacts_realized.json"
            self.assertFalse(
                realized.exists(),
                "realized sidecar must not be written for dry runs",
            )


if __name__ == "__main__":
    unittest.main()
