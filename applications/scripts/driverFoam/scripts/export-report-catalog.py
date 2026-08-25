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
#     export-report-catalog
#
# Description
#     Exports report metadata schemas to JSON.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Export the active plugin's report catalog to JSON.

Backend authors report definitions against ``openfoam_driver/core/report_catalog.py``'s
``ReportDefinition`` record; each plugin owns its own catalog (the built-in
cardiac plugin's lives at ``openfoam_driver/plugins/cardiacfoam/reports.py``)
and this script writes it to a stable JSON catalog for external consumers.
Defaults to the built-in cardiacFoam plugin, matching v1 behavior, unless
``--plugin`` selects otherwise.

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

from openfoam_driver.core.report_catalog import to_record  # noqa: E402


def build_catalog(plugin: str | None) -> dict:
    from openfoam_driver.core.plugin_interface import (  # noqa: E402
        default_driver_context,
        generic_openfoam_context,
        load_plugin_context,
    )

    if plugin == "none":
        context = generic_openfoam_context()
    elif plugin:
        context = load_plugin_context(plugin)
    else:
        context = default_driver_context()
    reports = context.capabilities.report_catalog.reports()
    return {
        "version": "1",
        "reports": [to_record(r) for r in reports],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True, help="output JSON path")
    parser.add_argument(
        "--plugin",
        help=(
            "Plugin whose report catalog to export: an installed plugin id, "
            "a trusted local-development import target "
            "(module.path:PluginClass), or 'none' for generic OpenFOAM. "
            "Defaults to built-in cardiacFoam."
        ),
    )
    args = parser.parse_args()
    catalog = build_catalog(args.plugin)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(catalog, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
