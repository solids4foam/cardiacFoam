#!/usr/bin/env python3
"""Check exported CPU-batched and CUDA-batched single-cell rates."""

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
    parser.add_argument("cpu_case", type=Path)
    parser.add_argument("cuda_case", type=Path)
    parser.add_argument(
        "--start-time", type=float, default=None,
        help="compare samples at or after this time (seconds)",
    )
    args = parser.parse_args()

    cpu_names, cpu = load_trace(args.cpu_case)
    cuda_names, cuda = load_trace(args.cuda_case)
    if args.start_time is not None:
        cpu = cpu[cpu[:, 0] >= args.start_time]
        cuda = cuda[cuda[:, 0] >= args.start_time]
        if not len(cpu) or not len(cuda):
            raise SystemExit("no samples at or after requested start time")
    if cpu_names != cuda_names or cpu.shape != cuda.shape:
        raise SystemExit("CPU and CUDA traces have different headers or shapes")
    names = [name for name in cpu_names if name.startswith("RATES_")]
    if not names:
        raise SystemExit("trace exports no RATES_* columns")
    indices = [cpu_names.index(name) for name in names]
    cpu_rates, cuda_rates = cpu[:, indices], cuda[:, indices]
    if not np.isfinite(cpu_rates).all() or not np.isfinite(cuda_rates).all():
        raise SystemExit("nonfinite rate values")
    cpu_nonzero = int(np.count_nonzero(cpu_rates))
    cuda_nonzero = int(np.count_nonzero(cuda_rates))
    max_abs = float(np.max(np.abs(cpu_rates-cuda_rates)))
    max_reference = float(np.max(np.abs(cpu_rates)))
    direct_names = [
        name for name in cpu_names
        if name not in {"time", "t"} and not name.startswith("AV_")
    ]
    direct_indices = [cpu_names.index(name) for name in direct_names]
    direct_diff = np.abs(cpu[:, direct_indices]-cuda[:, direct_indices])
    direct_max = float(np.max(direct_diff))
    print(f"samples={cpu.shape[0]} rate_columns={len(names)} "
          f"cpu_nonzero={cpu_nonzero} cuda_nonzero={cuda_nonzero} "
          f"max_abs_diff={max_abs:.9g} max_abs_rate={max_reference:.9g}")
    print(f"Vm/state/current/rate columns={len(direct_names)} "
          f"max_abs_diff_at_saved_precision={direct_max:.9g}")
    # The output format writes seven digits after the decimal point. This
    # tolerance allows the combined rounding of two backend traces while
    # still rejecting materially different exported derivatives.
    if not cpu_nonzero or not cuda_nonzero:
        print("Rates export: FAIL (one backend exported only zeros)")
        return 1
    if max_abs > 2.0e-7 or direct_max > 2.0e-7:
        print("Rates/backend parity: FAIL (difference exceeds output precision)")
        return 1
    print("Rates and direct state/current parity: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
