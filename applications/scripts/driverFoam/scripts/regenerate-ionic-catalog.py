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
#     regenerate-ionic-catalog
#
# Description
#     Regenerates state and algebraic properties in the ionic catalog.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Regenerate `states`, `algebraic`, `constants` tuples in ionic_model_catalog.py
from each model's `*_Names.H` header. All other catalogue fields (descriptions,
aliases, recommended_*, model_type, ...) are preserved untouched.

Idempotent: running this twice produces zero diff.

Usage:
    python applications/scripts/driverFoam/scripts/regenerate-ionic-catalog.py
    python applications/scripts/driverFoam/scripts/regenerate-ionic-catalog.py --check
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

DRIVER_PKG = Path(__file__).resolve().parents[1]
REPO_ROOT = DRIVER_PKG.parents[2]

# Make the openfoam_driver package importable when the script is run directly.
sys.path.insert(0, str(DRIVER_PKG))

from openfoam_driver.scripts._names_parser import (  # noqa: E402
    EXCLUDED_FROM_HEADER_SYNC,
    ParsedNames,
    find_names_header,
    parse_names_header,
)

CATALOG_PATH = (
    DRIVER_PKG
    / "openfoam_driver"
    / "plugins"
    / "cardiacfoam"
    / "ionic_model_catalog.py"
)
IONIC_MODELS_DIR = REPO_ROOT / "src" / "ionicModels"


def _format_tuple(items: tuple[str, ...]) -> str:
    if not items:
        return "()"
    return "(" + ", ".join(f'"{x}"' for x in items) + ")"


def rewrite_entry(text: str, model_name: str, parsed: ParsedNames) -> tuple[str, int]:
    """Replace the three data tuples for one model entry. Returns
    (new_text, n_substitutions). n_substitutions is 0 if the entry could not
    be matched (e.g. someone reformatted the catalogue).
    """
    pattern = re.compile(
        r'("'
        + re.escape(model_name)
        + r'":\s*IonicModelEntry\(\s*\n)'
        + r"(?P<i1>[ \t]+)states=\([^)]*\),[ \t]*\n"
        + r"(?P<i2>[ \t]+)algebraic=\([^)]*\),[ \t]*\n"
        + r"(?P<i3>[ \t]+)constants=\([^)]*\),"
    )

    def repl(m: re.Match[str]) -> str:
        head = m.group(1)
        i1, i2, i3 = m.group("i1"), m.group("i2"), m.group("i3")
        return (
            f"{head}"
            f"{i1}states={_format_tuple(parsed.states)},\n"
            f"{i2}algebraic={_format_tuple(parsed.algebraic)},\n"
            f"{i3}constants={_format_tuple(parsed.constants)},"
        )

    return pattern.subn(repl, text)


def collect_models() -> list[tuple[str, Path]]:
    """Return [(model_name, header_path), ...] for every model directory with
    a Names header, excluding the FDA manufactured semantic-label models.
    """
    out: list[tuple[str, Path]] = []
    for child in sorted(IONIC_MODELS_DIR.iterdir()):
        if not child.is_dir() or child.name == "lnInclude":
            continue
        if child.name in EXCLUDED_FROM_HEADER_SYNC:
            continue
        header = find_names_header(child)
        if header is None:
            continue
        out.append((child.name, header))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit non-zero if the catalogue would change; do not write.",
    )
    args = parser.parse_args()

    original = CATALOG_PATH.read_text(encoding="utf-8")
    text = original

    unmatched: list[str] = []
    for model_name, header in collect_models():
        parsed = parse_names_header(header)
        text, n = rewrite_entry(text, model_name, parsed)
        if n == 0:
            unmatched.append(model_name)
        elif n > 1:
            print(
                f"warning: {model_name} matched {n} times in catalogue "
                f"(expected 1)",
                file=sys.stderr,
            )

    if unmatched:
        print(
            "error: could not locate catalogue entries for:\n  "
            + "\n  ".join(unmatched)
            + "\n(catalogue file may have been reformatted - the regex expects "
            "states/algebraic/constants as the first three fields on single lines)",
            file=sys.stderr,
        )
        return 2

    if text == original:
        if not args.check:
            print("ionic_model_catalog.py: already in sync, no changes.")
        return 0

    if args.check:
        print(
            "ionic_model_catalog.py is OUT OF SYNC with *_Names.H headers.\n"
            "Run: python applications/scripts/driverFoam/scripts/"
            "regenerate-ionic-catalog.py",
            file=sys.stderr,
        )
        return 1

    CATALOG_PATH.write_text(text, encoding="utf-8")
    print(f"updated {CATALOG_PATH.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
