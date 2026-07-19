"""Tolerance diff of a fresh canonical convergence CSV against a committed reference.

Exit codes: 0 PASS, 1 FAIL, 2 SKIP (reference missing). Stdlib only.
Errors (L1/L2/Linf) compared by relative tolerance; rates by absolute tolerance.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

KEY = ("case", "variant", "dim", "N", "field")
ERR_COLS = ("L1", "L2", "Linf")
RATE_COLS = ("rate_L2", "rate_Linf")


def _load(path):
    with Path(path).open(newline="") as fh:
        return {tuple(r[k] for k in KEY): r for r in csv.DictReader(fh)}


def _rel_ok(a, b, rel):
    try:
        a, b = float(a), float(b)
    except (TypeError, ValueError):
        return (a or "") == (b or "")          # blank/non-numeric must match exactly
    denom = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / denom <= rel


def _abs_ok(a, b, tol):
    try:
        a, b = float(a), float(b)
    except (TypeError, ValueError):
        return (a or "") == (b or "")
    return abs(a - b) <= tol


def compare(fresh, reference, rel=0.05, rate_abs=0.15):
    """Return (ok, failures). Every reference row must have a matching fresh row
    within tolerance; extra fresh rows are ignored."""
    fresh_rows = _load(fresh)
    ref_rows = _load(reference)
    failures = []
    for key, ref in ref_rows.items():
        cur = fresh_rows.get(key)
        if cur is None:
            failures.append(f"missing row {key}")
            continue
        for col in ERR_COLS:
            if not _rel_ok(cur.get(col, ""), ref.get(col, ""), rel):
                failures.append(f"{key} {col}: fresh={cur.get(col)} ref={ref.get(col)} (rel>{rel})")
        for col in RATE_COLS:
            if not _abs_ok(cur.get(col, ""), ref.get(col, ""), rate_abs):
                failures.append(f"{key} {col}: fresh={cur.get(col)} ref={ref.get(col)} (abs>{rate_abs})")
    return (not failures), failures


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("fresh")
    ap.add_argument("reference")
    ap.add_argument("--rel", type=float, default=0.05)
    ap.add_argument("--rate-abs", type=float, default=0.15)
    args = ap.parse_args(argv)
    if not Path(args.reference).exists():
        print(f"SKIP: reference missing ({args.reference})")
        return 2
    ok, failures = compare(args.fresh, args.reference, args.rel, args.rate_abs)
    if ok:
        print(f"PASS: {args.fresh} matches {args.reference}")
        return 0
    print(f"FAIL: {len(failures)} mismatch(es) vs {args.reference}")
    for f in failures:
        print("  " + f)
    return 1


if __name__ == "__main__":
    sys.exit(main())
