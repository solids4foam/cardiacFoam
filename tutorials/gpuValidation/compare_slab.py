#!/usr/bin/env python3
"""Compare matched scalar and batched OpenFOAM slab output fields."""

import argparse
import re
from pathlib import Path

import numpy as np


FIELD_RE = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)",
    re.DOTALL,
)


def field(case, time, name, n_cells=None):
    path = case / time / name
    data = path.read_text()
    match = FIELD_RE.search(data)
    if match is None:
        uniform = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", data)
        if uniform is None or n_cells is None:
            raise ValueError(f"expected scalar internal field: {path}")
        return np.full(n_cells, float(uniform.group(1)))
    values = np.fromstring(match.group(2), sep=" ")
    if len(values) != int(match.group(1)) or not np.isfinite(values).all():
        raise ValueError(f"invalid or nonfinite field: {path}")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scalar_case", type=Path)
    parser.add_argument("batched_case", type=Path)
    parser.add_argument("--times", nargs="+", default=["0.005", "0.01", "0.015"])
    args = parser.parse_args()
    for time in args.times:
        scalar_vm = field(args.scalar_case, time, "Vm")
        batched_vm = field(args.batched_case, time, "Vm")
        scalar_activation = field(args.scalar_case, time, "activationTime", len(scalar_vm))
        batched_activation = field(args.batched_case, time, "activationTime", len(batched_vm))
        if scalar_vm.shape != batched_vm.shape or scalar_activation.shape != batched_activation.shape:
            raise ValueError(f"field sizes differ at {time}")
        vm_error_mv = 1000 * (batched_vm - scalar_vm)
        scalar_mask = scalar_activation > 0
        batched_mask = batched_activation > 0
        shared_mask = scalar_mask & batched_mask
        activation_difference_ms = (
            1000 * np.abs(batched_activation[shared_mask] - scalar_activation[shared_mask])
        )
        p95 = float(np.percentile(activation_difference_ms, 95)) if len(activation_difference_ms) else float("nan")
        print(
            f"{time}s cells={len(scalar_vm)} "
            f"activated_scalar={int(scalar_mask.sum())} "
            f"activated_batched={int(batched_mask.sum())} "
            f"Vm_RMSE_mV={np.sqrt(np.mean(vm_error_mv**2)):.6f} "
            f"Vm_max_mV={np.max(np.abs(vm_error_mv)):.6f} "
            f"activation_p95_ms={p95:.6f}"
        )


if __name__ == "__main__":
    main()
