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
#     test_utility_catalog_export
#
# Description
#     Tests utility catalog export logic and manifest coverage contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the utility-catalog exporter.

The utility catalog is an *agent tool catalog*: it is what lets a caller
build a valid invocation of an OpenFOAM utility and know which artifacts
come back. That makes ``positional_args``, ``flags`` (with their
``argument_kind``) and ``produces`` load-bearing, not decoration — an
earlier version of the exporter dropped all three, which is precisely why
the emitted JSON had no consumer.

Source of truth is the distributed ``utility.manifest.toml`` sidecar next to
each utility's C++, aggregated into ``core/utility_catalog.py``. These tests
pin the exporter to that dataclass so a new manifest field cannot be added
without also reaching the JSON.
"""

from __future__ import annotations

import dataclasses
import json
import subprocess
import sys
from pathlib import Path

from openfoam_driver.core.utility_catalog import UTILITY_CATALOG

REPO = Path(__file__).resolve().parents[3]
SCRIPT = REPO / "scripts" / "export-utility-catalog.py"


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


def test_exporter_writes_versioned_utility_list(tmp_path):
    data = _run(tmp_path / "u.json")
    assert data["version"] == "1"
    assert isinstance(data["utilities"], list)
    assert len(data["utilities"]) == len(UTILITY_CATALOG)


def test_every_catalog_entry_is_exported(tmp_path):
    data = _run(tmp_path / "u.json")
    assert {u["name"] for u in data["utilities"]} == set(UTILITY_CATALOG)


def test_entries_are_sorted_by_name(tmp_path):
    names = [u["name"] for u in _run(tmp_path / "u.json")["utilities"]]
    assert names == sorted(names), "output must be stable across runs"


# ---------------------------------------------------------------------------
# Nothing may be silently dropped
# ---------------------------------------------------------------------------


def test_no_manifest_field_is_dropped(tmp_path):
    """Every ``UtilityManifest`` field must reach the JSON.

    Guards the failure mode this exporter actually had: fields present in
    the TOML and in the dataclass, but absent from the exported record.
    """
    data = _run(tmp_path / "u.json")
    expected = {f.name for f in dataclasses.fields(next(iter(UTILITY_CATALOG.values())))}
    for record in data["utilities"]:
        assert set(record) == expected, (
            f"{record['name']}: exported fields do not match UtilityManifest; "
            f"missing={sorted(expected - set(record))} "
            f"extra={sorted(set(record) - expected)}"
        )


def test_flags_carry_the_full_flag_contract(tmp_path):
    data = _run(tmp_path / "u.json")
    flag_fields = {
        f.name
        for manifest in UTILITY_CATALOG.values()
        for flag in manifest.flags
        for f in dataclasses.fields(flag)
    }
    seen = False
    for record in data["utilities"]:
        for flag in record["flags"]:
            assert set(flag) == flag_fields
            seen = True
    assert seen, "no utility declares a flag; the assertion above never ran"


def test_produces_entries_survive_export(tmp_path):
    """``produces`` is how a caller knows what artifacts to expect."""
    data = _run(tmp_path / "u.json")
    exported = sum(len(u["produces"]) for u in data["utilities"])
    expected = sum(len(m.produces) for m in UTILITY_CATALOG.values())
    assert expected > 0, "catalog declares no produces; test is vacuous"
    assert exported == expected


def test_positional_args_survive_export(tmp_path):
    data = _run(tmp_path / "u.json")
    exported = sum(len(u["positional_args"]) for u in data["utilities"])
    expected = sum(len(m.positional_args) for m in UTILITY_CATALOG.values())
    assert expected > 0, "catalog declares no positional args; test is vacuous"
    assert exported == expected


# ---------------------------------------------------------------------------
# Portability
# ---------------------------------------------------------------------------


def test_source_paths_are_repo_relative(tmp_path):
    data = _run(tmp_path / "u.json")
    for record in data["utilities"]:
        assert not Path(record["source_path"]).is_absolute(), record["name"]


def test_output_has_no_python_repr_leaks(tmp_path):
    """No Python repr may reach the JSON (frozensets, Paths, tuples)."""
    raw = tmp_path / "u.json"
    _run(raw)
    text = raw.read_text()
    for leak in ("frozenset(", "PosixPath(", "WindowsPath("):
        assert leak not in text, f"{leak} leaked into the exported JSON"
