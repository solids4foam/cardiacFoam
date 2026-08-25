"""P2.9: end-to-end regression gate for the trust-boundary claims in
SECURITY.md.

One test per documented claim, exercised through the real entry points named in
that document -- the CLI ``run``/``step`` path and ``build_execution_inputs``
ingestion, not re-implemented logic. Claims SECURITY.md lists under
"Explicitly NOT mitigated" are covered too, as assertions that the *documented
current behaviour* is still accurate; the one open gap the document calls out
as a value-channel hole gets an explicit ``xfail(strict=True)`` asserting the
mitigated behaviour, so closing it forces this file (and SECURITY.md) to be
updated rather than leaving a stale claim behind.

Every test names the SECURITY.md line it gates in its docstring.
"""
from __future__ import annotations

import json
import os
import tempfile
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest import mock

import pytest

from openfoam_driver.tests.conftest import skip_without_monorepo

pytestmark = skip_without_monorepo

from openfoam_driver.cli import main
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.core.runtime.workflow_runner import (
    _resolve_case_cwd,
    _resolve_command,
    run_workflow_step,
)
from openfoam_driver.core.runtime.workflow_state import initial_workflow_state
from openfoam_driver.core.strict_planning import strict_plan

REPO_ROOT = Path(__file__).resolve().parents[6]
SECURITY_MD = REPO_ROOT / "applications" / "scripts" / "driverFoam" / "SECURITY.md"
SINGLE_CELL_ROOT = REPO_ROOT / "tutorials" / "electrophysiologyProtocols" / "singleCell"

CASE_NAME = "trustBoundaryCase"

# An Allrun that produces exactly the artifacts strict planning predicts for
# this case shape, so a permitted run reaches "completed" and a blocked one is
# distinguishable by the absence of these files.
ALLRUN_OK = (
    "#!/bin/sh\n"
    f"mkdir -p postProcessing 0.001\n"
    f"touch postProcessing/{CASE_NAME}_1.txt 0.001/Vm 0.001/AV_Ta\n"
    "exit 0\n"
)


def _write_case(root: Path, *, allrun: str = ALLRUN_OK, steps: list[dict] | None = None) -> Path:
    """Create a minimal runnable OpenFOAM case with an Allrun-owned workflow.

    Mirrors tests/core/test_cli_run_document.py::_write_case so these tests
    drive the same real planning path without needing a cardiacFoam binary.
    """
    case_root = root / CASE_NAME
    (case_root / "constant").mkdir(parents=True)
    (case_root / "system").mkdir()
    (case_root / "constant" / "electroProperties").write_text(
        (SINGLE_CELL_ROOT / "constant" / "electroProperties").read_text()
    )
    (case_root / "constant" / "physicsProperties").write_text(
        (SINGLE_CELL_ROOT / "constant" / "physicsProperties").read_text()
    )
    for name in ("controlDict", "fvSchemes", "fvSolution"):
        (case_root / "system" / name).write_text("\n")
    allrun_path = case_root / "Allrun"
    allrun_path.write_text(allrun)
    os.chmod(allrun_path, 0o755)
    del steps
    return case_root


def _spec_with_workflow(case_root: Path, *, steps: list[dict]) -> TutorialSpec:
    return TutorialSpec(
        name=case_root.name,
        case_root=case_root,
        setup_root=case_root,
        output_dir=case_root / "postProcessing",
        build_cases=lambda: [CaseConfig(case_id="default", params={})],
        apply_case=lambda *_args, **_kwargs: None,
        metadata={
            "entry_name": case_root.name,
            "entry_kind": "case_folder",
            "entry_path": case_root.name,
            "source_type": "filesystem_case",
            "workflow_family": None,
            "workflow_dag": {"steps": steps},
        },
    )


def _cli(argv: list[str]) -> tuple[int, dict]:
    """Invoke the real CLI and return (exit_code, parsed JSON report)."""
    out = StringIO()
    with redirect_stdout(out):
        code = main(argv)
    return code, json.loads(out.getvalue())


def _plan_to_file(tutorials_root: Path, doc_path: Path) -> dict:
    """`plan --strict --entry` and persist the emitted RunDocument."""
    code, report = _cli([
        "plan", "--strict", "--entry", CASE_NAME,
        "--tutorials-root", str(tutorials_root),
    ])
    assert code == 0, report
    run_document = report["run_document"]
    assert run_document is not None, report
    doc_path.write_text(json.dumps(run_document))
    return run_document


def _tampered_document(
    tutorials_root: Path,
    doc_path: Path,
    *,
    steps: list[dict] | None = None,
    launch: dict | None = None,
    config: dict | None = None,
) -> Path:
    """A real planned RunDocument with one field replaced by agent content.

    This is the adversarial shape SECURITY.md is written against: an agent
    hands back a document that is well-formed and passes every *other* gate,
    with exactly one hostile field. Starting from a genuinely planned document
    (rather than a hand-written stub) keeps the resulting diagnostics
    attributable to the tampered field alone, instead of drowning in unrelated
    config-validation errors.
    """
    document = _plan_to_file(tutorials_root, doc_path)
    if steps is not None:
        document["workflowDag"]["steps"] = steps
    if launch is not None:
        document["launch"] = {**document["launch"], **launch}
    if config is not None:
        document["config"] = config
    doc_path.write_text(json.dumps(document))
    return doc_path


def _hand_authored_document(
    doc_path: Path,
    *,
    case_root: Path,
    steps: list[dict],
    launch: dict | None = None,
    config: dict | None = None,
) -> Path:
    """A RunDocument built by hand, never derived from `plan --strict`.

    Complements `_tampered_document` (a real planned document with one field
    swapped): this constructs the whole document from scratch -- the shape
    an agent could hand back without ever invoking the planner. Mirrors the
    inline hand-authored documents in
    tests/core/test_cli_run_document.py::test_run_document_rejects_unknown_command.

    `config` defaults to a real, validation-passing config harvested by
    planning a throwaway case of the same tutorial family in an isolated
    directory, so tests targeting the workflow/launch boundary aren't also
    incidentally blocked by unrelated physics-config validation errors.
    """
    if config is None:
        with tempfile.TemporaryDirectory() as seed_dir:
            seed_root = Path(seed_dir)
            _write_case(seed_root)
            seed_document = _plan_to_file(seed_root, seed_root / "seed.json")
            config = seed_document["config"]
    document = {
        "version": "3",
        "id": "hand-authored",
        "name": "hand-authored",
        "status": "planned",
        "config": config,
        "launch": launch if launch is not None else {
            "caseRoot": str(case_root),
            "outputDir": str(case_root / "driverfoam-output"),
        },
        "workflowDag": {
            "schema_version": "1",
            "step_status_values": [
                "pending", "running", "completed", "failed", "skipped",
            ],
            "steps": steps,
        },
    }
    doc_path.write_text(json.dumps(document))
    return doc_path


def _step(command: str, **overrides) -> dict:
    step = {
        "id": "s", "command": command, "args": [], "cwd": ".",
        "depends_on": [], "produces": [], "consumes": [],
        "retry_policy": {}, "command_display": command,
    }
    step.update(overrides)
    return step


def _codes(payload: dict) -> set[str]:
    return {d.get("code") for d in payload.get("diagnostics", [])}


def _plan_codes(report: dict) -> set[str]:
    """Every diagnostic code in a `plan --strict` report, across its buckets."""
    codes: set[str] = set()
    for key, value in report.items():
        if key.endswith("_diagnostics") and isinstance(value, list):
            codes.update(d.get("code") for d in value if isinstance(d, dict))
    return codes


# --------------------------------------------------------------------------
# Documented-closed claim: "launch.caseRoot ... when DRIVERFOAM_ALLOWED_RUNS_ROOT
# is set, both must resolve under it" (Trust boundaries / Mitigations).
# --------------------------------------------------------------------------

def test_case_root_outside_allowed_runs_root_is_rejected_before_execution() -> None:
    """SECURITY.md: "opt-in DRIVERFOAM_ALLOWED_RUNS_ROOT containment"."""
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(tutorials_root)
        doc_path = tutorials_root / "run.json"
        _plan_to_file(tutorials_root, doc_path)

        elsewhere = tutorials_root / "allowed-elsewhere"
        elsewhere.mkdir()
        with mock.patch.dict(
            os.environ, {"DRIVERFOAM_ALLOWED_RUNS_ROOT": str(elsewhere)}
        ):
            code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code != 0, payload
        assert "case_root_outside_allowed_root" in _codes(payload), payload
        # Rejected at ingestion: the case script never ran.
        assert not (case_root / "0.001").exists()


def test_output_dir_outside_allowed_runs_root_is_rejected_before_execution() -> None:
    """SECURITY.md: containment applies to `outputDir`, not just `caseRoot`."""
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        allowed = root / "allowed"
        allowed.mkdir()
        case_root = _write_case(allowed)
        doc_path = root / "run.json"
        _plan_to_file(allowed, doc_path)

        document = json.loads(doc_path.read_text())
        document["launch"]["outputDir"] = str(root / "outside-results")
        doc_path.write_text(json.dumps(document))

        with mock.patch.dict(
            os.environ, {"DRIVERFOAM_ALLOWED_RUNS_ROOT": str(allowed)}
        ):
            code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code != 0, payload
        assert "output_dir_outside_allowed_root" in _codes(payload), payload
        assert not (case_root / "0.001").exists()


def test_symlinked_case_root_cannot_escape_allowed_runs_root() -> None:
    """SECURITY.md: paths are "resolved to canonical absolute paths", so
    containment "runs on resolved paths, so a symlink ... cannot escape"."""
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        allowed = root / "allowed"
        allowed.mkdir()
        outside = root / "outside"
        outside.mkdir()
        real_case = _write_case(outside)
        doc_path = root / "run.json"
        _plan_to_file(outside, doc_path)

        # A symlink that *lexically* sits inside the allowed root but points out.
        link = allowed / CASE_NAME
        link.symlink_to(real_case)
        document = json.loads(doc_path.read_text())
        document["launch"]["caseRoot"] = str(link)
        document["launch"]["outputDir"] = str(link / "out")
        doc_path.write_text(json.dumps(document))

        with mock.patch.dict(
            os.environ, {"DRIVERFOAM_ALLOWED_RUNS_ROOT": str(allowed)}
        ):
            code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code != 0, payload
        assert "case_root_outside_allowed_root" in _codes(payload), payload
        assert not (real_case / "0.001").exists()


def test_allowed_runs_root_permits_a_contained_case() -> None:
    """Containment is a boundary, not a blanket refusal: an in-root case runs.

    Without this, every rejection test above would also pass against a driver
    that refused everything.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(tutorials_root)
        doc_path = tutorials_root / "run.json"
        _plan_to_file(tutorials_root, doc_path)

        with mock.patch.dict(
            os.environ, {"DRIVERFOAM_ALLOWED_RUNS_ROOT": str(tutorials_root)}
        ):
            code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code == 0, payload
        assert payload["status"] == "ok"
        assert (case_root / "0.001" / "Vm").exists()


# --------------------------------------------------------------------------
# Documented-closed claim: "launch.caseRoot must be an existing, runnable
# OpenFOAM case (registry._case_is_runnable)".
# --------------------------------------------------------------------------

def test_case_root_must_be_an_existing_runnable_openfoam_case() -> None:
    """SECURITY.md: "caseRoot must be a runnable OpenFOAM case"."""
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(tutorials_root)

        not_a_case = tutorials_root / "just-a-directory"
        not_a_case.mkdir()
        doc_path = _tampered_document(
            tutorials_root, tutorials_root / "run.json",
            launch={"caseRoot": str(not_a_case), "outputDir": str(not_a_case / "out")},
        )
        code, payload = _cli(["run", "--run-document", str(doc_path)])
        assert code != 0, payload
        assert "case_root_not_a_runnable_case" in _codes(payload), payload

        missing = tutorials_root / "does-not-exist"
        doc_path = _tampered_document(
            tutorials_root, tutorials_root / "run2.json",
            launch={"caseRoot": str(missing), "outputDir": str(missing / "out")},
        )
        code, payload = _cli(["run", "--run-document", str(doc_path)])
        assert code != 0, payload
        assert "case_root_missing" in _codes(payload), payload


# --------------------------------------------------------------------------
# Documented-closed claim: the command allowlist (validate_workflow_commands).
# --------------------------------------------------------------------------

def test_command_authorization_rejects_an_unauthorized_bare_command() -> None:
    """SECURITY.md: only "a known OpenFOAM/driver command, an Allrun-family
    case script, a registered utility, or an installed OpenFOAM app"."""
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(tutorials_root)
        sentinel = tutorials_root / "PWNED"
        doc_path = _tampered_document(
            tutorials_root, tutorials_root / "run.json",
            steps=[_step("touch", args=[str(sentinel)])],
        )
        code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code != 0, payload
        assert "unknown_workflow_command" in _codes(payload), payload
        assert not sentinel.exists(), "unauthorized command must not execute"
        assert not (case_root / "0.001").exists()


def test_command_authorization_rejects_an_absolute_path_command() -> None:
    """SECURITY.md: "Absolute-path ... commands are rejected"."""
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(tutorials_root)
        sentinel = tutorials_root / "PWNED"
        doc_path = _tampered_document(
            tutorials_root, tutorials_root / "run.json",
            steps=[_step("/usr/bin/touch", args=[str(sentinel)])],
        )
        code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code != 0, payload
        assert "unknown_workflow_command" in _codes(payload), payload
        message = " ".join(
            d["message"] for d in payload["diagnostics"]
            if d.get("code") == "unknown_workflow_command"
        )
        assert "explicit path" in message, payload
        assert not sentinel.exists()


def test_command_authorization_rejects_an_arbitrary_relative_script() -> None:
    """SECURITY.md: "arbitrary ./script commands are rejected" -- only
    ./Allrun-family case scripts may be given in path form."""
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(tutorials_root)
        sentinel = tutorials_root / "PWNED"
        pwn = case_root / "pwn.sh"
        pwn.write_text(f"#!/bin/sh\ntouch {sentinel}\n")
        os.chmod(pwn, 0o755)

        doc_path = _tampered_document(
            tutorials_root, tutorials_root / "run.json",
            steps=[_step("./pwn.sh")],
        )
        code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code != 0, payload
        assert "unknown_workflow_command" in _codes(payload), payload
        assert not sentinel.exists()

        # ...while the documented exception (./Allrun-family) is accepted.
        doc_path = _tampered_document(
            tutorials_root, tutorials_root / "run_allrun.json",
            steps=[_step("./Allrun", id="run", produces=[])],
        )
        code, payload = _cli(["run", "--run-document", str(doc_path)])
        assert code == 0, payload
        assert (case_root / "0.001" / "Vm").exists()


def test_command_allowlist_has_one_owner_shared_by_both_producers() -> None:
    """SECURITY.md: "Single command allowlist owner
    (validate_workflow_commands), enforced once at ingestion" -- so the strict
    planner and the run-document adapter cannot drift apart.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        # Producer 1: strict planning from a driver-owned Python workflow.
        planner_root = Path(temp_dir) / "planner"
        planner_root.mkdir()
        planner_case = _write_case(planner_root)
        with mock.patch(
            "openfoam_driver.core.strict_planning.load_entry_spec",
            return_value=_spec_with_workflow(
                planner_case,
                steps=[{"id": "run", "command": "curl", "depends_on": []}],
            ),
        ):
            plan_report = strict_plan(
                CASE_NAME,
                overrides={"tutorials_root": str(planner_root)},
            ).to_json()
        assert "unknown_workflow_command" in _plan_codes(plan_report), plan_report

        # Producer 2: an agent-authored RunDocument carrying the same command,
        # ingested through the run path instead of the planner.
        adapter_root = Path(temp_dir) / "adapter"
        adapter_root.mkdir()
        _write_case(adapter_root)
        doc_path = _tampered_document(
            adapter_root, adapter_root / "run.json", steps=[_step("curl")],
        )
        run_code, run_payload = _cli(["run", "--run-document", str(doc_path)])
        assert run_code != 0, run_payload
        assert "unknown_workflow_command" in _codes(run_payload), run_payload


# --------------------------------------------------------------------------
# Documented-closed claim: "No command shadowing: bare names resolve via PATH
# only; only Allrun-family resolve case-locally".
# --------------------------------------------------------------------------

def test_case_directory_cannot_shadow_a_trusted_path_binary() -> None:
    """SECURITY.md: "No command shadowing"."""
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(
            tutorials_root,
            steps=[{"id": "mesh", "command": "blockMesh", "depends_on": []}],
        )
        sentinel = tutorials_root / "PWNED"
        shadow = case_root / "blockMesh"
        shadow.write_text(f"#!/bin/sh\ntouch {sentinel}\n")
        os.chmod(shadow, 0o755)

        # The resolver is the enforcement point: a case-local `blockMesh` is
        # never picked up, while an Allrun-family name deliberately is.
        assert _resolve_command("blockMesh", case_root) == "blockMesh"
        assert _resolve_command("Allrun", case_root) == str(case_root / "Allrun")

        # End-to-end: running the step never executes the case-local shadow,
        # whether or not a real blockMesh exists on this machine's PATH.
        doc_path = _hand_authored_document(
            tutorials_root / "run.json",
            case_root=case_root,
            steps=[_step("blockMesh", id="mesh")],
        )
        _cli(["run", "--run-document", str(doc_path)])
        assert not sentinel.exists(), "case-local binary shadowed a PATH command"


# --------------------------------------------------------------------------
# Documented-closed claim: "Workflow cwd cannot escape caseRoot".
# --------------------------------------------------------------------------

def test_workflow_cwd_cannot_escape_case_root() -> None:
    """SECURITY.md: "Workflow `cwd` cannot escape `caseRoot`".

    Two layers: ingestion refuses a lexically-escaping `cwd`, and the runner
    re-checks the *resolved* path so a symlink inside the case cannot escape
    either.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(tutorials_root)

        # Layer 1 -- ingestion, via the real CLI run path.
        doc_path = _hand_authored_document(
            tutorials_root / "run.json",
            case_root=case_root,
            steps=[_step("Allrun", cwd="../..")],
        )
        code, payload = _cli(["run", "--run-document", str(doc_path)])
        assert code != 0, payload
        assert "workflow_cwd_not_case_relative" in _codes(payload), payload

        # Layer 2 -- the runner's own resolved-path check.
        with pytest.raises(ValueError, match="escapes case root"):
            _resolve_case_cwd(case_root, "../..")

        outside = tutorials_root / "outside"
        outside.mkdir()
        (case_root / "escape").symlink_to(outside)
        with pytest.raises(ValueError, match="escapes case root"):
            _resolve_case_cwd(case_root, "escape")


# --------------------------------------------------------------------------
# Documented-closed claim: "Steps run argv-style (no shell)".
# --------------------------------------------------------------------------

def test_steps_run_argv_style_so_arguments_are_not_shell_interpreted() -> None:
    """SECURITY.md: "Steps run argv-style (no shell)"."""
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        sentinel = tutorials_root / "PWNED"
        recorded = "args.txt"
        case_root = _write_case(
            tutorials_root,
            allrun=(
                "#!/bin/sh\n"
                f"printf '%s' \"$1\" > {recorded}\n"
                f"mkdir -p postProcessing 0.001\n"
                f"touch postProcessing/{CASE_NAME}_1.txt 0.001/Vm 0.001/AV_Ta\n"
                "exit 0\n"
            ),
        )
        doc_path = _hand_authored_document(
            tutorials_root / "run.json",
            case_root=case_root,
            steps=[_step("Allrun", args=[f"; touch {sentinel}"])],
        )

        code, payload = _cli(["run", "--run-document", str(doc_path)])
        assert code == 0, payload
        assert not sentinel.exists(), "argument was interpreted by a shell"
        assert (case_root / recorded).read_text() == f"; touch {sentinel}"


# --------------------------------------------------------------------------
# Documented-closed claim: "config via validate_run" at ingestion.
# --------------------------------------------------------------------------

def test_invalid_config_blocks_execution_at_ingestion() -> None:
    """SECURITY.md: the RunDocument's `config` is validated at ingestion
    ("`config` via `validate_run`") before anything is executed.

    The malformed value here is an out-of-enum string (a domain-semantic
    violation `validate_run` is meant to catch), not a wrong Python type for
    a phase slice. The wrong-type case is gated separately by
    `test_non_mapping_config_phase_blocks_execution_at_ingestion` below;
    keeping the two apart preserves this test's original claim (domain
    semantics) as its own regression gate.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(tutorials_root)
        doc_path = _hand_authored_document(
            tutorials_root / "run.json",
            case_root=case_root,
            steps=[_step("Allrun")],
            config={
                "anatomy": {},
                "physics": {"type": "monodomain", "myocardiumSolver": "not-a-real-solver"},
                "stimulus": {},
                "solver": {},
            },
        )
        code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code != 0, payload
        assert any(
            d.get("level") == "error" for d in payload.get("diagnostics", ())
        ), payload
        assert not (case_root / "0.001").exists(), "config was not gated before execution"


def test_non_mapping_config_phase_blocks_execution_at_ingestion() -> None:
    """A wrong-*type* config phase must produce a diagnostic, not a crash.

    P2.2 opened the core schema's `config` to `additionalProperties: true`
    with no per-phase type constraint, so a phase value can legally be any
    JSON type by the time it reaches `validate_run`. Before the guard in
    `specs/validation.py::_non_mapping_phase_errors`, a non-dict phase
    (`{"anatomy": "not-an-object"}`) reached `_flatten_context` and raised an
    uncaught `AttributeError` through the real
    `driverFoam run --run-document` path -- a traceback instead of the
    diagnostic SECURITY.md promises. This is the regression gate for that
    fix: the CLI must exit non-zero with a parseable JSON payload.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(tutorials_root)
        doc_path = _hand_authored_document(
            tutorials_root / "run.json",
            case_root=case_root,
            steps=[_step("Allrun")],
            config={
                "anatomy": "not-an-object",
                "physics": {},
                "stimulus": {},
                "solver": {},
            },
        )
        # A traceback escaping `main` would fail here before any assertion.
        code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code != 0, payload
        diagnostics = payload.get("diagnostics", ())
        assert any(d.get("level") == "error" for d in diagnostics), payload
        assert any(
            d.get("field") == "anatomy" and "must be an object" in d.get("message", "")
            for d in diagnostics
        ), payload
        assert not (case_root / "0.001").exists(), "config was not gated before execution"


# --------------------------------------------------------------------------
# Documented-OPEN claims. These assert that SECURITY.md's "Explicitly NOT
# mitigated" section is an accurate statement about today's code, so that
# closing one of these holes shows up here as a failure and forces the
# document to be updated instead of silently going stale.
# --------------------------------------------------------------------------

def test_case_script_caveat_is_still_documented_and_still_true() -> None:
    """SECURITY.md documents case scripts as "untrusted, unsandboxed by
    design" -- "running a case runs its code". Regression-check that this is
    still an accurate statement about current behavior, not a stale claim.
    """
    text = SECURITY_MD.read_text()
    assert "unsandboxed by design" in text
    assert "Arbitrary code inside an invoked `Allrun`" in text

    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        # The Allrun writes *outside* its own case root: nothing confines it.
        escaped = tutorials_root / "written-by-allrun.txt"
        _write_case(
            tutorials_root,
            allrun=(
                "#!/bin/sh\n"
                f"printf 'allrun ran unsandboxed' > {escaped}\n"
                f"mkdir -p postProcessing 0.001\n"
                f"touch postProcessing/{CASE_NAME}_1.txt 0.001/Vm 0.001/AV_Ta\n"
                "exit 0\n"
            ),
        )
        doc_path = tutorials_root / "run.json"
        _plan_to_file(tutorials_root, doc_path)
        code, payload = _cli(["run", "--run-document", str(doc_path)])

        assert code == 0, payload
        assert escaped.read_text() == "allrun ran unsandboxed", (
            "SECURITY.md still claims Allrun contents are unsandboxed; if this "
            "assertion fails a sandbox was added -- update SECURITY.md."
        )


def test_run_workflow_step_is_still_a_trusted_unvalidating_primitive() -> None:
    """SECURITY.md: "run_workflow_step is a trusted low-level primitive: a
    Python caller that invokes it directly with an unvalidated case_root /
    ... / command bypasses path and command validation."

    Asserts the documented caveat is still true. If validation is ever added to
    the runner, this fails and SECURITY.md must be corrected.
    """
    text = SECURITY_MD.read_text()
    assert "run_workflow_step` is a trusted low-level primitive" in text

    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        # Not an OpenFOAM case at all, and a command no allowlist authorizes.
        bare_dir = root / "not-a-case"
        bare_dir.mkdir()
        sentinel = root / "written-by-unvalidated-runner.txt"
        dag = {
            "schema_version": "1",
            "step_status_values": [
                "pending", "running", "completed", "failed", "skipped",
            ],
            "steps": [{
                "id": "s", "command": "touch", "args": [str(sentinel)],
                "cwd": ".", "depends_on": [], "produces": [], "consumes": [],
                "retry_policy": {}, "command_display": "touch",
            }],
        }
        result = run_workflow_step(
            dag, initial_workflow_state(dag), "s",
            case_root=bare_dir, log_dir=root / "logs",
        )
        assert result.state.steps[0].status == "completed"
        assert sentinel.exists(), (
            "SECURITY.md still claims run_workflow_step performs no path/command "
            "validation; if this fails, validation was added -- update SECURITY.md."
        )


def test_override_values_containing_a_coded_entry_are_rejected() -> None:
    """Asserts the *mitigated* behaviour: `step --apply` refuses an override
    whose value smuggles executable OpenFOAM code into a case dictionary.

    `mutators._format_value` (tier 1, the path almost every override takes)
    now rejects any override value containing `#`, `;`, or a newline before
    it is ever written to a case dictionary file. See SECURITY.md.

    Real entry point: `driverFoam step --run-document <doc> --step <id> --apply
    <overrides.json>`, which routes through specs.apply_overrides.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(tutorials_root)
        control_dict = case_root / "system" / "controlDict"
        control_dict.write_text("deltaT    0.001;\nendTime    1;\n")
        doc_path = tutorials_root / "run.json"
        _plan_to_file(tutorials_root, doc_path)

        payload = '#codeStream { code #{ os << 0.001; #}; }'
        overrides = tutorials_root / "ov.json"
        overrides.write_text(json.dumps([{"driver_path": "deltaT", "value": payload}]))

        code, report = _cli([
            "step", "--run-document", str(doc_path), "--step", "run",
            "--apply", str(overrides),
        ])

        assert code != 0, report
        assert "#codeStream" not in control_dict.read_text()
