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
import sys

import check_against_reference as _checker

# Reuse the numeric checker's own row loader so the key is EXACTLY its
# KEY = ("case","variant","dim","N","field") — including N. References carry
# multiple N-rows per (case,variant,dim,field) group, so N must stay in the
# key or an extra/missing resolution row would slip through — the very gap
# this gate exists to close. Calling _load (not re-listing columns) makes
# drift from the checker impossible.
KEY = _checker.KEY


def keyset(csv_path) -> set[tuple]:
    return set(_checker._load(csv_path).keys())


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
