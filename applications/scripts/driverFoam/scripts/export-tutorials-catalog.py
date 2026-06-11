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
#     export-tutorials-catalog
#
# Description
#     Exports tutorial metadata definitions to JSON format.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Export the tutorial catalog.

Cross-checks ``tutorials_display.TUTORIALS`` against
``core.runtime.registry.REGISTERED_TUTORIALS`` so we cannot ship a
tutorial catalog entry without a backend factory, or omit a registered
tutorial. Writes a stable JSON shape to the requested output path.

Source-of-truth stays in Python (the user confirmed: "the tutorials are
currently running based on hardcoded scripts in the backend. that is okay").
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openfoam_driver.core.runtime.registry import REGISTERED_TUTORIALS  # noqa: E402
from openfoam_driver.tutorials_display import TUTORIALS, to_record  # noqa: E402


def build_catalog() -> dict:
    display_ids = {t.id for t in TUTORIALS}
    registry_ids = set(REGISTERED_TUTORIALS)

    only_in_display = display_ids - registry_ids
    only_in_registry = registry_ids - display_ids
    if only_in_display or only_in_registry:
        raise SystemExit(
            "tutorials_display.TUTORIALS is out of sync with "
            "REGISTERED_TUTORIALS. "
            f"only-in-display={sorted(only_in_display)} "
            f"only-in-registry={sorted(only_in_registry)}. "
            "Either add a TutorialDisplay row or remove it; both "
            "sets must match exactly."
        )

    # Stable order: keep REGISTERED_TUTORIALS' declared order so the
    # JSON is reproducible regardless of TUTORIALS tuple order.
    by_id = {t.id: t for t in TUTORIALS}
    return {
        "version": "1",
        "tutorials": [to_record(by_id[name]) for name in REGISTERED_TUTORIALS],
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
