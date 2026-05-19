#!/usr/bin/env python3
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

from openfoam_driver.utility_catalog import UTILITY_CATALOG  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[4]


def _manifest_to_record(manifest) -> dict:
    """Serialise a ``UtilityManifest`` to a JSON-ready dict."""
    return {
        "name": manifest.name,
        "description": manifest.description,
        "purpose": manifest.purpose,
        "inputs": list(manifest.inputs),
        "outputs": list(manifest.outputs),
        "requires_mesh": manifest.requires_mesh,
        "flags": [
            {
                "name": f.name,
                "description": f.description,
                "takes_value": f.takes_value,
            }
            for f in manifest.flags
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
