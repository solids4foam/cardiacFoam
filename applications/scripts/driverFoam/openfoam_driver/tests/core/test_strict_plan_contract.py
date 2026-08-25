"""``plan --strict`` must always answer in its own contract.

Every failure mode driverFOAM knows about is reported as a JSON document with
``status: "failed"`` and structured diagnostics -- that document IS the
interface an agent consumes. A failure that escapes as an unhandled exception
gives the caller a traceback on stderr and nothing on stdout, so
``json.load()`` raises and the agent has to parse English prose to find out
what went wrong.

``myocardiumSolver`` is the selector the whole coeffs scope name is derived
from, so a case missing it genuinely cannot be planned. Failing is correct;
failing outside the contract is not.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

DRIVER_ROOT = Path(__file__).resolve().parents[3]
_SINGLE_CELL = (
    DRIVER_ROOT.parents[2] / "tutorials" / "electrophysiologyProtocols" / "singleCell"
)


def _case_without_solver_selector(tmp_path: Path) -> Path:
    tutorials_root = tmp_path / "tutorials"
    case = tutorials_root / "case"
    shutil.copytree(_SINGLE_CELL, case)
    ep = case / "constant" / "electroProperties"
    ep.write_text(
        ep.read_text().replace(
            "myocardiumSolver singleCellSolver;", "myocardiumSolvr singleCellSolver;"
        )
    )
    return tutorials_root


def test_unresolvable_solver_selector_is_reported_as_a_diagnostic(tmp_path):
    from openfoam_driver.core.strict_planning import strict_plan

    tutorials_root = _case_without_solver_selector(tmp_path)
    payload = strict_plan(
        "case",
        entry_kind="case_folder",
        overrides={"tutorials_root": str(tutorials_root)},
        openfoam_bashrc="/no/such/openfoam/bashrc",
    ).to_json()

    assert payload["status"] == "failed"
    # Assert the contract, not the bucket: what matters is that the reason is
    # machine-readable somewhere in the document. "missing_solver" is the code
    # the plugin already uses for this, so no new vocabulary is invented.
    reported = [
        item
        for key, value in payload.items()
        if key.endswith("_diagnostics") or key == "catalog_coverage_errors"
        for item in (value or [])
        if item.get("code") == "missing_solver"
    ]
    assert reported, f"no missing_solver diagnostic in {list(payload)}"
    assert any("myocardiumSolver" in item["message"] for item in reported)


def test_the_cli_still_emits_parseable_json(tmp_path):
    """The end-to-end shape an agent actually consumes."""
    tutorials_root = _case_without_solver_selector(tmp_path)
    result = subprocess.run(
        [
            sys.executable, "-m", "openfoam_driver", "plan", "--strict",
            "--tutorials-root", str(tutorials_root),
            "--entry", "case", "--entry-kind", "case_folder",
            "--openfoam-bashrc", "/no/such/openfoam/bashrc",
        ],
        cwd=DRIVER_ROOT, capture_output=True, text=True,
    )
    assert result.returncode != 0, "an unplannable case must still exit non-zero"
    payload = json.loads(result.stdout)  # must not raise
    assert payload["status"] == "failed"
