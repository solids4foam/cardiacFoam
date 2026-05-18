"""Tests for the tutorials-catalog exporter.

The user confirmed tutorials live as hardcoded Python in the backend
("the tutorials are currently running based on hardcoded scripts in
the backend. that is okay. We can pass them to json if that is not
done yet"). The hardcoded list is
``openfoam_driver.core.runtime.registry.REGISTERED_TUTORIALS``; this
exporter serializes it (plus a thin display-metadata layer) to
JSON for external consumers.

Source-of-truth stays in Python.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
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
    from openfoam_driver.core.runtime.registry import REGISTERED_TUTORIALS

    data = _run(tmp_path / "t.json")
    exported = {t["id"] for t in data["tutorials"]}
    registered = set(REGISTERED_TUTORIALS)
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
