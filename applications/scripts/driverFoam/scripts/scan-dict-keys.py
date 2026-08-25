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
#     scan-dict-keys
#
# Description
#     Generates drift reports for dictionary keys and parameter schemas.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Dict-keys drift report for the cardiacFoam driverFoam Python catalogue.

Compares dictionary keys actually read by the C++ source tree against
``driver_path`` entries declared in ``dict_entries.py``.

Accuracy is ~80%; false positives/negatives are acceptable — this is an
ad-hoc human-review tool, not a CI contract test.

Usage:
    python applications/scripts/driverFoam/scripts/scan-dict-keys.py
    python applications/scripts/driverFoam/scripts/scan-dict-keys.py --limit 10

Exit code is always 0.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

DRIVER_PKG = Path(__file__).resolve().parents[1]
REPO_ROOT = DRIVER_PKG.parents[2]

# Make the openfoam_driver package importable when the script is run directly.
sys.path.insert(0, str(DRIVER_PKG))

from openfoam_driver.scripts._dict_keys_scanner import (  # noqa: E402
    CataloguePath,
    DictRead,
    iter_catalogue_paths,
    scan_dict_reads,
    strict_dict_key_report,
)

def _plugin_scan_inputs(plugin: str | None):
    """Resolve the active plugin's C++ mapping and dictionary catalogue.

    Mirrors ``core/strict_planning.py``: the source roots, the allowlist and
    the catalogue are plugin-owned provenance, not constants of this script.
    Returns ``(cxx_mapping | None, entries)``.
    """
    from openfoam_driver.core.plugin_interface import (
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
    mapping = context.capabilities.cxx_mapping.profile().cxx_mapping
    entries = tuple(context.capabilities.dictionaries.entries())
    return mapping, entries


# OpenFOAM FoamFile-header boilerplate keys that are never user-facing.
IGNORED_KEYS: frozenset[str] = frozenset(
    {
        "version",
        "format",
        "class",
        "object",
        "location",
        "dimensions",
        "internalField",
        "boundaryField",
        "FoamFile",
        "note",
        "arch",
        "root",
        "case",
        "time",
        "path",
    }
)


# ---------------------------------------------------------------------------
# Helpers


def _short_path(p: Path) -> str:
    """Return path relative to repo root for display."""
    try:
        return str(p.relative_to(REPO_ROOT))
    except ValueError:
        return str(p)


def _site_label(read: DictRead) -> str:
    return f"{_short_path(read.source_file)}:{read.line}"


# ---------------------------------------------------------------------------
# Report sections


def _section_a(
    code_keys: dict[str, list[DictRead]],
    cat_leaves: set[str],
    limit: int,
) -> tuple[int, str]:
    """Keys read in C++ but absent from the catalogue (leaf comparison)."""
    absent = {k: v for k, v in code_keys.items() if k not in cat_leaves and k not in IGNORED_KEYS}

    if not absent:
        body = "  (none — all C++ keys matched a catalogue leaf name)\n"
        return 0, body

    col_w = max(len(k) for k in absent) + 2

    lines: list[str] = []
    lines.append(
        f"  {'key':<{col_w}}| count | example sites"
    )
    lines.append(
        f"  {'-' * col_w}+-------+{'-' * 45}"
    )

    shown = 0
    for key, reads in sorted(absent.items(), key=lambda kv: -len(kv[1])):
        if shown >= limit:
            lines.append(f"  ... ({len(absent) - shown} more rows — use --limit to raise cap)")
            break
        sites = ", ".join(_site_label(r) for r in reads[:3])
        if len(reads) > 3:
            sites += f", ... (+{len(reads) - 3} more)"
        lines.append(f"  {key:<{col_w}}|  {len(reads):>4} | {sites}")
        shown += 1

    return len(absent), "\n".join(lines) + "\n"


def _section_b(
    cat_paths: list[CataloguePath],
    code_keys_set: set[str],
    limit: int,
) -> tuple[int, str]:
    """Catalogue driver_paths with no apparent C++ reader (leaf comparison)."""
    stale: list[CataloguePath] = []
    for cp in cat_paths:
        # Skip entries whose leaf is a wildcard placeholder — they cannot appear
        # as literal strings in code.
        if cp.has_wildcard and cp.dynamic_path:
            continue
        if cp.leaf not in code_keys_set:
            stale.append(cp)

    if not stale:
        body = "  (none — every catalogue leaf name was found in C++ reads)\n"
        return 0, body

    col_w = max(len(cp.driver_path) for cp in stale) + 2

    lines: list[str] = []
    lines.append(f"  {'driver_path':<{col_w}}| note")
    lines.append(f"  {'-' * col_w}+{'-' * 50}")

    shown = 0
    for cp in stale:
        if shown >= limit:
            lines.append(f"  ... ({len(stale) - shown} more rows — use --limit to raise cap)")
            break
        note = f"leaf '{cp.leaf}' not found in any dict read"
        if cp.has_wildcard:
            note += " (has wildcard segments)"
        lines.append(f"  {cp.driver_path:<{col_w}}| {note}")
        shown += 1

    return len(stale), "\n".join(lines) + "\n"


def _section_c(
    code_subdicts: dict[str, list[DictRead]],
    cat_parent_segs: set[str],
    limit: int,
) -> tuple[int, str]:
    """Sub-dict names found in C++ vs. parent-path segments in the catalogue."""
    all_names = sorted(
        set(code_subdicts.keys()) | cat_parent_segs
    )

    if not all_names:
        body = "  (no subDict calls or catalogue parent segments found)\n"
        return 0, body

    col_w = max(len(n) for n in all_names) + 2

    lines: list[str] = []
    lines.append(
        f"  {'subdict name':<{col_w}}| in code  | in catalogue"
    )
    lines.append(
        f"  {'-' * col_w}+----------+--------------"
    )

    unmatched = 0
    shown = 0
    for name in all_names:
        if shown >= limit:
            remaining = len(all_names) - shown
            lines.append(f"  ... ({remaining} more rows — use --limit to raise cap)")
            break
        in_code = "✓" if name in code_subdicts else "✗"
        in_cat = "✓" if name in cat_parent_segs else "✗"
        if in_code != in_cat:
            unmatched += 1
        lines.append(f"  {name:<{col_w}}|    {in_code}     |    {in_cat}")
        shown += 1

    # Count remaining unmatched beyond the cap.
    for name in all_names[shown:]:
        in_code = name in code_subdicts
        in_cat = name in cat_parent_segs
        if in_code != in_cat:
            unmatched += 1

    return unmatched, "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=50,
        metavar="N",
        help="Maximum rows per section (default: 50).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Fail on drift not covered by the reviewed allowlist, and fail "
            "when allowlist entries become unused. Note that "
            "'driverFoam plan --strict' already runs this same check and "
            "reports it as plugin_dict_key_* diagnostics."
        ),
    )
    parser.add_argument(
        "--plugin",
        help=(
            "Plugin whose C++ mapping and catalogue to scan: an installed "
            "plugin id, a trusted local-development import target "
            "(module.path:PluginClass), or 'none' for generic OpenFOAM. "
            "Defaults to built-in cardiacFoam."
        ),
    )
    args = parser.parse_args()

    mapping, entries = _plugin_scan_inputs(args.plugin)
    if mapping is None:
        print(
            "Active plugin declares no C++ mapping; nothing to scan.",
            file=sys.stderr,
        )
        return 0

    src_roots = [root for root in mapping.source_roots if root.is_dir()]
    for root in mapping.source_roots:
        if not root.is_dir():
            print(f"Source root unavailable, skipped: {root}", file=sys.stderr)
    if not src_roots:
        print("No usable C++ source root; nothing to scan.", file=sys.stderr)
        return 0

    if args.strict:
        status = 0
        for root in src_roots:
            report = strict_dict_key_report(
                root,
                allowlist_path=mapping.allowlist_path,
                entries=entries,
            )
            print(json.dumps(report.to_json(), indent=2))
            if report.status != "ok":
                status = 1
        return status

    for root in src_roots:
        print(f"Scanning {_short_path(root)} ...", file=sys.stderr)
    reads = [read for root in src_roots for read in scan_dict_reads(root)]
    cat_paths = list(iter_catalogue_paths(entries))

    # Index C++ reads.
    key_reads: dict[str, list[DictRead]] = defaultdict(list)
    subdict_reads: dict[str, list[DictRead]] = defaultdict(list)
    for r in reads:
        if r.kind == "key":
            key_reads[r.name].append(r)
        else:
            subdict_reads[r.name].append(r)

    code_keys_set: set[str] = set(key_reads.keys())
    code_subdicts: set[str] = set(subdict_reads.keys())  # type: ignore[assignment]

    # Index catalogue.
    cat_leaves: set[str] = set()
    cat_parent_segs: set[str] = set()
    for cp in cat_paths:
        if not (cp.has_wildcard and cp.dynamic_path):
            cat_leaves.add(cp.leaf)
        for seg in cp.parents:
            if not _WILDCARD_RE_SIMPLE.match(seg):
                cat_parent_segs.add(seg)

    # Build sections.
    count_a, body_a = _section_a(dict(key_reads), cat_leaves, args.limit)
    count_b, body_b = _section_b(cat_paths, code_keys_set, args.limit)
    count_c, body_c = _section_c(dict(subdict_reads), cat_parent_segs, args.limit)

    # Print report.
    print()
    print("== Section A: dictionary keys read in src/ but absent from catalogue ==")
    print("(leaf-name comparison; ~80% accuracy — manually verify before catalogue edits)")
    print()
    print(body_a)

    print("== Section B: catalogue driver_paths with no apparent C++ reader ==")
    print("(leaf-name comparison; may be false positives if the key is read indirectly)")
    print()
    print(body_b)

    print("== Section C: sub-dict names ==")
    print("(comparison between subDict(\"X\") calls and catalogue parent paths)")
    print()
    print(body_c)

    print("== Summary ==")
    print(f"  C++ unique keys read:                   {len(code_keys_set)}")
    print(f"  Catalogue unique leaf names:            {len(cat_leaves)}")
    print(f"  Section A count (potential gaps):       {count_a}")
    print(f"  Section B count (potential stale):      {count_b}")
    print(f"  Section C unmatched:                    {count_c}")
    print()

    return 0


# Simple wildcard segment detector (no re import needed at module level,
# but we need it here).
import re as _re  # noqa: E402
_WILDCARD_RE_SIMPLE = _re.compile(r"<[^>]+>")


if __name__ == "__main__":
    sys.exit(main())
