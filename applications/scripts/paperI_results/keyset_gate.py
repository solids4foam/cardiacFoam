"""Exact key-set gate for a fresh canonical convergence CSV vs. a reference.

Complements check_against_reference.py's numeric tolerance diff, which only
requires every *reference* row to have a matching fresh row (extra fresh rows
are silently ignored). This gate additionally requires the fresh and
reference key sets to be exactly equal, so stale/unexpected rows in the fresh
CSV cannot coexist with a PASS.

Exit codes: 0 key sets equal, 1 key sets differ. Stdlib only.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import check_against_reference as _checker

# Group on the same columns check_against_reference.py keys rows by, minus
# "N": this gate checks which (case, variant, dim, field) combinations exist
# at all, independent of mesh resolution. Deriving from _checker.KEY (rather
# than re-listing the column names) keeps the shared columns from drifting
# out of sync with the numeric checker.
KEY = tuple(k for k in _checker.KEY if k != "N")


def keyset(csv_path) -> set[tuple]:
    with Path(csv_path).open(newline="") as fh:
        return {tuple(r[k] for k in KEY) for r in csv.DictReader(fh)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fresh")
    ap.add_argument("reference")
    args = ap.parse_args(argv)

    fresh_keys = keyset(args.fresh)
    ref_keys = keyset(args.reference)

    if fresh_keys == ref_keys:
        print(f"PASS: key sets match ({len(fresh_keys)} keys, {KEY})")
        return 0

    only_fresh = fresh_keys - ref_keys
    only_ref = ref_keys - fresh_keys
    print(f"FAIL: key sets differ ({KEY})")
    for k in sorted(only_fresh):
        print(f"  only in fresh:     {k}")
    for k in sorted(only_ref):
        print(f"  only in reference: {k}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
