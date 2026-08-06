#!/usr/bin/env python3
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
# Script
#     export-dict-catalog
#
# Description
#     Exports dictionary parameters and catalog items to JSON.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Export ``dict_entries`` and the ionic / active-tension catalogs to JSON.

The exporter reads ``DictEntry`` objects and fans each one out into one record
per phase the entry declares. An entry tagged ``phases={"anatomy", "physics"}``
therefore appears in BOTH the ``anatomy`` and ``physics`` buckets; each emitted
record carries a single ``phase`` field equal to its bucket and preserves the
full ``phases`` list for validation and agent consumers.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openfoam_driver.dict_entries import (  # noqa: E402
    PHYSICS_PROPERTY_ENTRIES,
    get_electro_property_entry_groups,
)
from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG  # noqa: E402
from openfoam_driver.plugins.cardiacfoam.active_tension_catalog import (  # noqa: E402
    ACTIVE_TENSION_MODEL_CATALOG,
)

PHASES = ("anatomy", "physics", "stimulus", "solver")


def _all_entries():
    """Yield every ``DictEntry`` known to the backend, regardless of group."""
    entries: list[DictEntry] = list(PHYSICS_PROPERTY_ENTRIES)
    for group in get_electro_property_entry_groups().values():
        entries.extend(group)
    yield from entries


def _entry_to_record(e) -> dict:
    """Flatten a ``DictEntry`` to a JSON record.

    ``phases`` arrives as a ``frozenset`` (unordered, not JSON-serialisable);
    we emit a sorted list so the catalog JSON is stable across runs.
    """
    d = asdict(e)
    d["phases"] = sorted(e.phases)
    return d


def build_catalog() -> dict:
    by_phase: dict[str, list] = {p: [] for p in PHASES}
    for e in _all_entries():
        if not e.phases:
            raise SystemExit(f"entry missing phases: {e.driver_path}")
        record = _entry_to_record(e)
        # Fan-out: one record per declared phase. Each emitted record is
        # stamped with a single `phase` (its bucket) while keeping `phases`
        # so downstream consumers know the entry's other homes.
        for ph in e.phases:
            if ph not in by_phase:
                raise SystemExit(
                    f"entry {e.driver_path} has unknown phase {ph!r}"
                )
            by_phase[ph].append({**record, "phase": ph})
    return {
        "version": "1",
        "phases": {p: {"entries": by_phase[p]} for p in PHASES},
        "ionic_models": [asdict(m) for m in IONIC_MODEL_CATALOG.values()],
        "active_tension_models": [
            asdict(m) for m in ACTIVE_TENSION_MODEL_CATALOG.values()
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True, help="output JSON path")
    args = parser.parse_args()
    catalog = build_catalog()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(catalog, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
