"""Canonical Paper I convergence-results schema and helpers.

Single source of truth for the CSV layout every conforming case writes to
<case>/setup/results/<case>_convergence.csv. Reconcile these columns with the
results output added to post_processing_manufactured.py when those edits land.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

CANONICAL_FIELDS = [
    "case", "variant", "dim", "N", "h",
    "field", "L1", "L2", "Linf", "rate_L2", "rate_Linf",
]


def observed_order(e_prev, e_cur, h_prev, h_cur):
    """p = log(e_prev/e_cur) / log(h_prev/h_cur); None on degenerate input."""
    try:
        e_prev, e_cur = float(e_prev), float(e_cur)
        h_prev, h_cur = float(h_prev), float(h_cur)
    except (TypeError, ValueError):
        return None
    if min(e_prev, e_cur, h_prev, h_cur) <= 0 or h_prev == h_cur:
        return None
    return math.log(e_prev / e_cur) / math.log(h_prev / h_cur)


def _group_key(row):
    return (row["case"], row["variant"], row["dim"], row["field"])


def fill_rates(rows):
    """Return rows (copied) with rate_L2/rate_Linf filled per (case,variant,dim,
    field) group, coarse->fine by descending h; blank on the coarsest member."""
    out = [dict(r) for r in rows]
    groups: dict = {}
    for r in out:
        groups.setdefault(_group_key(r), []).append(r)
    for members in groups.values():
        members.sort(key=lambda r: float(r["h"]), reverse=True)
        prev = None
        for r in members:
            if prev is None:
                r["rate_L2"] = r["rate_Linf"] = ""
            else:
                p2 = observed_order(prev.get("L2"), r.get("L2"), prev["h"], r["h"])
                pi = observed_order(prev.get("Linf"), r.get("Linf"), prev["h"], r["h"])
                r["rate_L2"] = "" if p2 is None else f"{p2:.2f}"
                r["rate_Linf"] = "" if pi is None else f"{pi:.2f}"
            prev = r
    out.sort(key=lambda r: (r["case"], r["variant"], r["dim"],
                            r["field"], -float(r["h"])))
    return out


def write_canonical(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CANONICAL_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in CANONICAL_FIELDS})
