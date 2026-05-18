#!/usr/bin/env python3
"""Export ``report_catalog.REPORTS`` to JSON.

Backend authors report definitions in ``openfoam_driver/report_catalog.py``;
this script writes them to a stable JSON catalog for external consumers.

URL templates are emitted verbatim — substitution of ``{port}`` and
``{kind}`` happens outside this exporter. The Python side never knows the
runtime port, which keeps 4Dpapers swappable.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openfoam_driver.report_catalog import REPORTS, to_record  # noqa: E402


def build_catalog() -> dict:
    return {
        "version": "1",
        "reports": [to_record(r) for r in REPORTS],
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
