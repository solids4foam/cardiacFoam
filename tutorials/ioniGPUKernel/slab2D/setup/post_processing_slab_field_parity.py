#!/usr/bin/env python3
"""Report voltage and activation parity for two matching slab outputs.

The first case is the reference; the second is the candidate. Use this for
scalar-versus-batched or host-batched-versus-CUDA-batched field comparisons.
"""

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
        values = np.full(n_cells, float(uniform.group(1)))
        if not np.isfinite(values).all():
            raise ValueError(f"nonfinite field: {path}")
        return values
    values = np.fromstring(match.group(2), sep=" ")
    if len(values) != int(match.group(1)) or not np.isfinite(values).all():
        raise ValueError(f"invalid or nonfinite field: {path}")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_case", type=Path)
    parser.add_argument("candidate_case", type=Path)
    parser.add_argument(
        "--times", nargs="+", default=["0.005", "0.01", "0.015"],
        help="written time directories to compare, in seconds",
    )
    args = parser.parse_args()
    for time in args.times:
        reference_vm = field(args.reference_case, time, "Vm")
        candidate_vm = field(args.candidate_case, time, "Vm")
        reference_activation = field(
            args.reference_case, time, "activationTime", len(reference_vm)
        )
        candidate_activation = field(
            args.candidate_case, time, "activationTime", len(candidate_vm)
        )
        if (
            reference_vm.shape != candidate_vm.shape
            or reference_activation.shape != candidate_activation.shape
        ):
            raise ValueError(f"field sizes differ at {time}")
        vm_error_mv = 1000 * (candidate_vm - reference_vm)
        reference_mask = reference_activation > 0
        candidate_mask = candidate_activation > 0
        shared_mask = reference_mask & candidate_mask
        activation_difference_ms = (
            1000 * np.abs(
                candidate_activation[shared_mask] - reference_activation[shared_mask]
            )
        )
        p95 = (
            float(np.percentile(activation_difference_ms, 95))
            if len(activation_difference_ms)
            else float("nan")
        )
        print(
            f"{time}s cells={len(reference_vm)} "
            f"activated_scalar={int(reference_mask.sum())} "
            f"activated_batched={int(candidate_mask.sum())} "
            f"activation_mask_mismatch={int(np.count_nonzero(reference_mask != candidate_mask))} "
            f"Vm_RMSE_mV={np.sqrt(np.mean(vm_error_mv**2)):.6f} "
            f"Vm_max_mV={np.max(np.abs(vm_error_mv)):.6f} "
            f"activation_p95_ms={p95:.6f}"
        )


if __name__ == "__main__":
    main()
