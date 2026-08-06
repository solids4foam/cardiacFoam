"""Compare two one-row interface-metric CSV files."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_row(path: Path) -> dict[str, str]:
    with path.open() as handle:
        return next(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-10)
    parser.add_argument("--relative-tolerance", type=float, default=1e-8)
    args = parser.parse_args()

    reference = read_row(args.reference)
    candidate = read_row(args.candidate)
    if reference.keys() != candidate.keys():
        missing = sorted(reference.keys() - candidate.keys())
        extra = sorted(candidate.keys() - reference.keys())
        raise SystemExit(f"CSV schema mismatch: missing={missing}, extra={extra}")
    rows = []
    failures = []
    for key, reference_text in reference.items():
        if key in {"method", "assembly", "fieldSource", "time"}:
            continue
        reference_value = float(reference_text)
        candidate_value = float(candidate[key])
        difference = abs(candidate_value - reference_value)
        tolerance = args.absolute_tolerance + args.relative_tolerance*abs(reference_value)
        passed = math.isfinite(difference) and difference <= tolerance
        rows.append((key, reference_value, candidate_value, difference, tolerance, passed))
        if not passed:
            failures.append(key)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("metric", "serial", "parallel", "absoluteDifference", "tolerance", "pass"))
        writer.writerows(rows)

    if failures:
        raise SystemExit("metric mismatch: " + ", ".join(failures))
    print(f"PASS: {len(rows)} metrics agree")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
