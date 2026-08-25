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
#     export-utility-catalog
#
# Description
#     Exports utility metadata to JSON.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Export the utility manifest catalog to JSON.

Each entry in the catalog is serialised to a flat JSON record; the ``flags``
list is inlined as an array of objects. ``source_path`` is converted to a
string relative to the repository root so the output is portable.

Usage::

    python scripts/export-utility-catalog.py --out /tmp/utility-catalog.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openfoam_driver.core.utility_catalog import UTILITY_CATALOG  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[4]


def _manifest_to_record(manifest) -> dict:
    """Serialise a ``UtilityManifest`` to a JSON-ready dict.

    Every field of ``UtilityManifest`` is emitted. The catalog is an agent
    tool-catalog: ``positional_args`` and ``produces`` are what let a caller
    build an invocation and know what artifacts come back, so dropping them
    (as an earlier version did) left the JSON unable to serve that purpose.
    """
    return {
        "name": manifest.name,
        "description": manifest.description,
        "purpose": manifest.purpose,
        "inputs": list(manifest.inputs),
        "outputs": list(manifest.outputs),
        "requires_mesh": manifest.requires_mesh,
        "positional_args": [
            {
                "name": a.name,
                "argument_kind": a.argument_kind,
                "description": a.description,
            }
            for a in manifest.positional_args
        ],
        "flags": [
            {
                "name": f.name,
                "description": f.description,
                "takes_value": f.takes_value,
                "argument_kind": f.argument_kind,
                "required": f.required,
                "default": f.default,
            }
            for f in manifest.flags
        ],
        "produces": [
            {
                "artifact_id": pr.artifact_id,
                "path_pattern": pr.path_pattern,
                "format": pr.format,
                "description": pr.description,
                "produced_by": pr.produced_by,
                "variables": list(pr.variables),
                "optional": pr.optional,
                "time_indexed": pr.time_indexed,
            }
            for pr in manifest.produces
        ],
        "example": manifest.example,
        "category": manifest.category,
        "source_path": str(
            manifest.source_path.relative_to(REPO_ROOT)
        ),
    }


def build_catalog() -> dict:
    return {
        "version": "1",
        "utilities": [
            _manifest_to_record(m)
            for m in sorted(UTILITY_CATALOG.values(), key=lambda m: m.name)
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
    print(f"Wrote {len(catalog['utilities'])} utility entries to {out_path}")


if __name__ == "__main__":
    main()
