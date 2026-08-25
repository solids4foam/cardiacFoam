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
#     test_report_catalog_export
#
# Description
#     Tests report catalog export logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the report-catalog exporter.

Backend Python declares report *definitions* (id, title, kind, URL template,
the ``applicable_when`` predicate, and a default-visible flag); the exporter
serializes them to JSON. The Python side never substitutes ``{port}`` or
``{kind}``, so 4Dpapers (or any future report backend) is swappable without
re-exporting.

For v1 the URL template is ``http://localhost:{port}/{kind}`` (no
``{runId}`` — the user's 4Dpapers app does not route by run yet). The
``applicable_when`` predicate language is flat key-equality:
``None`` ⇒ always applicable, ``{"phase.field": value}`` ⇒ AND of
equality checks. Anything else is invalid and the v1 evaluator throws.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SCRIPT = REPO / "scripts" / "export-report-catalog.py"


def _run(out_path: Path) -> dict:
    subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(out_path)],
        cwd=REPO,
        check=True,
    )
    return json.loads(out_path.read_text())


# ---------------------------------------------------------------------------
# Schema / shape
# ---------------------------------------------------------------------------


def test_exporter_writes_versioned_report_list(tmp_path):
    data = _run(tmp_path / "r.json")
    assert data["version"] == "1"
    assert isinstance(data["reports"], list)
    assert len(data["reports"]) >= 1, "v1 must ship at least the stub entry"
    sample = data["reports"][0]
    expected = {
        "id",
        "title",
        "kind",
        "url_template",
        "applicable_when",
        "show_by_default",
    }
    assert expected <= sample.keys(), (
        f"missing keys: {expected - sample.keys()}"
    )


def test_stub_entry_present(tmp_path):
    """v1 must ship the offline-bundled stub entry so the section
    works even when 4Dpapers is not running."""
    data = _run(tmp_path / "r.json")
    ids = {r["id"] for r in data["reports"]}
    assert "stub" in ids, f"expected 'stub' entry; got {sorted(ids)}"
    stub = next(r for r in data["reports"] if r["id"] == "stub")
    assert stub["url_template"] == "/reports/stub.html"
    assert stub["kind"] == "iframe"


def test_v1_url_templates_have_no_runid_placeholder(tmp_path):
    """4Dpapers v1 does not route by run; runId is a v2 concern.
    Any leaked ``{runId}`` here would force consumers to invent one."""
    data = _run(tmp_path / "r.json")
    for r in data["reports"]:
        assert "{runId}" not in r["url_template"], (
            f"report {r['id']!r} leaked {{runId}} into v1 url_template"
        )


def test_remote_entries_use_4dpapers_template(tmp_path):
    """At least one non-stub entry should target the 4Dpapers backend
    so the section is actually wired to the real app, not just the
    bundled fallback."""
    data = _run(tmp_path / "r.json")
    remote = [r for r in data["reports"] if r["id"] != "stub"]
    assert remote, "expected at least one 4Dpapers-backed report entry"
    for r in remote:
        assert r["url_template"].startswith("http://localhost:{port}/"), (
            f"report {r['id']!r} url_template does not target the "
            f"4Dpapers backend: {r['url_template']!r}"
        )

# ---------------------------------------------------------------------------
# applicable_when evaluator (the v1 predicate language)
# ---------------------------------------------------------------------------


def test_applicable_when_none_matches_anything():
    from openfoam_driver.core.report_catalog import matches

    assert matches(None, {"physics": {"ionic_model": "TenTusscher"}}) is True


def test_applicable_when_flat_equality_matches():
    from openfoam_driver.core.report_catalog import matches

    pred = {"physics.ionic_model": "TenTusscher"}
    cfg = {"physics": {"ionic_model": "TenTusscher"}}
    assert matches(pred, cfg) is True


def test_applicable_when_flat_equality_rejects_mismatch():
    from openfoam_driver.core.report_catalog import matches

    pred = {"physics.ionic_model": "TenTusscher"}
    cfg = {"physics": {"ionic_model": "FentonKarma"}}
    assert matches(pred, cfg) is False


def test_applicable_when_multi_key_is_AND():
    from openfoam_driver.core.report_catalog import matches

    pred = {
        "physics.ionic_model": "TenTusscher",
        "anatomy.mesh": "biventricular",
    }
    cfg = {
        "physics": {"ionic_model": "TenTusscher"},
        "anatomy": {"mesh": "biventricular"},
    }
    assert matches(pred, cfg) is True
    cfg["anatomy"]["mesh"] = "single-cell"
    assert matches(pred, cfg) is False


def test_applicable_when_missing_path_is_not_a_match():
    from openfoam_driver.core.report_catalog import matches

    pred = {"physics.ionic_model": "TenTusscher"}
    cfg = {"physics": {}}
    assert matches(pred, cfg) is False


def test_applicable_when_unknown_operator_raises():
    """v2 may add operators; v1 must refuse silently-mis-filtering."""
    import pytest

    from openfoam_driver.core.report_catalog import matches

    pred = {"physics.ionic_model": {"$in": ["TenTusscher"]}}
    with pytest.raises(ValueError, match="unsupported"):
        matches(pred, {"physics": {"ionic_model": "TenTusscher"}})
