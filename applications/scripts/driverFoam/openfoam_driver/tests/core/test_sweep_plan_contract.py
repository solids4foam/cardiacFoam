"""``sweep-plan`` must answer in its own contract, like ``plan --strict``.

A sweep already reports a per-case failure structurally -- ``status:
"failed"`` plus ``materialization_error`` on that case -- and keeps going.
That is the whole point of a sweep: one bad axis value should cost you one
case, not the command.

Two paths bypassed it. The per-case guard caught only ``OSError`` and
``ValueError``, so anything a tutorial factory raised (a ``KeyError`` for an
unknown ionic model, say) escaped and took the run down with a traceback --
zero bytes on stdout for the caller. And the spec file was read before any
guard at all, so a malformed JSON spec did the same.
"""

from __future__ import annotations

import json
from pathlib import Path

from openfoam_driver.core.runtime.sweep_runner import sweep_plan

_SPEC = {
    "base": {
        "entry": "singleCell",
        "case_dir_name": "electrophysiologyProtocols/singleCell",
        "setup_dir_name": "setup",
    },
    "sweep": {
        "mode": "zip",
        "independent": {
            "ionic_model": ["Courtemanche", "TotallyFakeModel"],
            "tissue": ["myocyte", "myocyte"],
        },
        "dependent": [
            {"name": "caseId", "derive": "case_id_template",
             "of": ["ionic_model", "tissue"]},
            {"name": "output_dir_name", "derive": "output_dir_name_template",
             "of": ["ionic_model", "tissue"]},
        ],
    },
}


def test_a_factory_failure_fails_one_case_not_the_command(tmp_path):
    """An unknown ionic model must cost exactly one case."""
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text(json.dumps(_SPEC))

    report = sweep_plan(spec_path, output_dir=tmp_path / "out")

    failed = [c for c in report["cases"] if c["status"] == "failed"]
    assert len(failed) == 1, report["cases"]
    assert "TotallyFakeModel" in failed[0]["materialization_error"]
    # ...and the good case still planned.
    assert any(c["status"] != "failed" for c in report["cases"]), report["cases"]


def test_a_malformed_spec_is_reported_structurally(tmp_path):
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text("{ not json")

    report = sweep_plan(spec_path, output_dir=tmp_path / "out")

    assert report["case_count"] == 0
    assert report["cases"] == []
    assert "spec_error" in report
    json.dumps(report)  # must be serialisable


def test_a_valid_spec_reports_no_spec_error(tmp_path):
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text(json.dumps(_SPEC))
    report = sweep_plan(spec_path, output_dir=tmp_path / "out")
    assert "spec_error" not in report


def test_a_malformed_spec_still_exits_non_zero(tmp_path):
    """Structured is not the same as successful."""
    import subprocess
    import sys

    spec_path = tmp_path / "sweep.json"
    spec_path.write_text("{ not json")
    driver_root = Path(__file__).resolve().parents[3]

    result = subprocess.run(
        [
            sys.executable, "-m", "openfoam_driver", "sweep-plan",
            "--spec", str(spec_path), "--output-dir", str(tmp_path / "out"),
        ],
        cwd=driver_root, capture_output=True, text=True,
    )
    assert result.returncode != 0, result.stdout
    assert json.loads(result.stdout)["spec_error"]
