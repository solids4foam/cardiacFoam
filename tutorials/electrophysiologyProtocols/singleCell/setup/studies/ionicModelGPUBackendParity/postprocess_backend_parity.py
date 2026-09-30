#!/usr/bin/env python3
"""Check exported rates, states, and currents from two batched backends."""

import argparse
from pathlib import Path

import numpy as np


def load_trace(case: Path) -> tuple[list[str], np.ndarray]:
    files = list((case / "postProcessing").glob("*.txt"))
    if len(files) != 1:
        raise SystemExit(f"expected one trace in {case}/postProcessing, found {len(files)}")
    with files[0].open() as stream:
        header = stream.readline().split()
    values = np.loadtxt(files[0], skiprows=1)
    if values.ndim == 1:
        values = values[None, :]
    if values.shape[1] != len(header):
        raise SystemExit(f"header/data width mismatch in {files[0]}")
    return header, values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_case", type=Path)
    parser.add_argument("candidate_case", type=Path)
    parser.add_argument(
        "--start-time", type=float, default=None,
        help="compare samples at or after this time (seconds)",
    )
    args = parser.parse_args()

    reference_names, reference = load_trace(args.reference_case)
    candidate_names, candidate = load_trace(args.candidate_case)
    if args.start_time is not None:
        reference = reference[reference[:, 0] >= args.start_time]
        candidate = candidate[candidate[:, 0] >= args.start_time]
        if not len(reference) or not len(candidate):
            raise SystemExit("no samples at or after requested start time")
    if reference_names != candidate_names or reference.shape != candidate.shape:
        raise SystemExit("reference and candidate traces have different headers or shapes")
    names = [name for name in reference_names if name.startswith("RATES_")]
    if not names:
        raise SystemExit("trace exports no RATES_* columns")
    indices = [reference_names.index(name) for name in names]
    reference_rates, candidate_rates = reference[:, indices], candidate[:, indices]
    if not np.isfinite(reference_rates).all() or not np.isfinite(candidate_rates).all():
        raise SystemExit("nonfinite rate values")
    reference_nonzero = int(np.count_nonzero(reference_rates))
    candidate_nonzero = int(np.count_nonzero(candidate_rates))
    max_abs = float(np.max(np.abs(reference_rates-candidate_rates)))
    max_reference = float(np.max(np.abs(reference_rates)))
    direct_names = [
        name for name in reference_names
        if name not in {"time", "t"} and not name.startswith("AV_")
    ]
    direct_indices = [reference_names.index(name) for name in direct_names]
    direct_diff = np.abs(reference[:, direct_indices]-candidate[:, direct_indices])
    direct_max = float(np.max(direct_diff))
    print(f"samples={reference.shape[0]} rate_columns={len(names)} "
          f"reference_nonzero={reference_nonzero} candidate_nonzero={candidate_nonzero} "
          f"max_abs_diff={max_abs:.9g} max_abs_rate={max_reference:.9g}")
    print(f"Vm/state/current/rate columns={len(direct_names)} "
          f"max_abs_diff_at_saved_precision={direct_max:.9g}")
    # The output format writes seven digits after the decimal point. This
    # tolerance allows the combined rounding of two backend traces while
    # still rejecting materially different exported derivatives.
    if not reference_nonzero or not candidate_nonzero:
        print("Rates export: FAIL (one backend exported only zeros)")
        return 1
    if max_abs > 2.0e-7 or direct_max > 2.0e-7:
        print("Rates/backend parity: FAIL (difference exceeds output precision)")
        return 1
    print("Rates and direct state/current backend parity: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
