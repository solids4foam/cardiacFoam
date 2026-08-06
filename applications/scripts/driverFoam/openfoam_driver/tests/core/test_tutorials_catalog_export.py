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
#     test_tutorials_catalog_export
#
# Description
#     Tests tutorials catalog export logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the tutorials-catalog exporter.

The user confirmed tutorials live as hardcoded Python in the backend
("the tutorials are currently running based on hardcoded scripts in
the backend. that is okay. We can pass them to json if that is not
done yet"). The hardcoded list is
``openfoam_driver.core.runtime.registry.list_tutorials()``; this
exporter serializes it (plus a thin display-metadata layer) to
JSON for external consumers.

Source-of-truth stays in Python.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SCRIPT = REPO / "scripts" / "export-tutorials-catalog.py"


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


def test_exporter_writes_versioned_tutorial_list(tmp_path):
    data = _run(tmp_path / "t.json")
    assert data["version"] == "1"
    assert isinstance(data["tutorials"], list)
    assert len(data["tutorials"]) >= 1
    sample = data["tutorials"][0]
    expected = {"id", "title", "summary", "thumbnail", "tags", "preset"}
    assert expected <= sample.keys(), (
        f"missing keys: {expected - sample.keys()}"
    )


def test_every_registered_tutorial_is_exported(tmp_path):
    """The display catalog cannot ship a tutorial card whose backend
    factory does not exist, and cannot omit a registered tutorial.
    Either drift makes the catalog lie about what the backend can run."""
    from openfoam_driver.core.runtime.registry import list_tutorials, _normalized_registry

    data = _run(tmp_path / "t.json")
    exported = {t["id"] for t in data["tutorials"]}
    registered = set(list_tutorials())
    assert exported == registered, (
        f"exported vs registered mismatch — only-in-exported: "
        f"{exported - registered}, only-in-registered: "
        f"{registered - exported}"
    )


def test_preset_uses_flat_key_equality(tmp_path):
    """The ``preset`` shape mirrors the v1 ``applicableWhen`` predicate:
    flat ``\"phase.field\": value`` keys. Nested dicts (operator-style)
    are reserved for v2 and must not leak into v1 output."""
    data = _run(tmp_path / "t.json")
    for t in data["tutorials"]:
        for key, value in t["preset"].items():
            assert "." in key, (
                f"tutorial {t['id']!r} preset key {key!r} is not "
                f"dotted phase.field form"
            )
            assert not isinstance(value, dict), (
                f"tutorial {t['id']!r} preset {key!r} uses an operator "
                f"object {value!r}; v1 is flat key-equality only"
            )


def test_tags_are_strings(tmp_path):
    data = _run(tmp_path / "t.json")
    for t in data["tutorials"]:
        assert isinstance(t["tags"], list)
        for tag in t["tags"]:
            assert isinstance(tag, str)
