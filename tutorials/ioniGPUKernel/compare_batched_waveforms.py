#!/usr/bin/env python3
"""Compare CPU/CUDA batched Vm traces with activation-time alignment."""

import argparse
import math
from pathlib import Path

import numpy as np

from slab2D.setup.post_processing_slab_field_parity import field


ACTIVATION_P95_LIMIT_MS = 0.1
ACTIVATED_COUNT_LIMIT_FRACTION = 0.01
ACTIVATED_COUNT_LIMIT_FLOOR = 5
PEAK_P95_LIMIT_MV = 2.0
APD90_ABS_LIMIT_MS = 2.0
APD90_REL_LIMIT = 0.02
ALIGNED_OUTSIDE_UPSTROKE_RMSE_LIMIT_MV = 0.5
UPSTROKE_HALF_WIDTH_S = 0.002


def times_in(case: Path) -> dict[float, str]:
    found = {}
    for child in case.iterdir():
        if not child.is_dir():
            continue
        try:
            numeric = float(child.name)
        except ValueError:
            continue
        if (child / "Vm").is_file():
            found[numeric] = child.name
    return found


def apd90(times: np.ndarray, traces: np.ndarray, activation: np.ndarray,
          active: np.ndarray) -> np.ndarray:
    result = np.full(traces.shape[1], np.nan, dtype=float)
    for cell in np.flatnonzero(active):
        start = max(0, int(np.searchsorted(times, activation[cell], side="right") - 1))
        tail = traces[start:, cell]
        if tail.size < 3:
            continue
        peak_rel = int(np.argmax(tail))
        peak_idx = start + peak_rel
        peak = traces[peak_idx, cell]
        rest = np.min(traces[:peak_idx + 1, cell])
        if peak <= rest:
            continue
        threshold = rest + 0.1*(peak - rest)
        crossings = np.flatnonzero(traces[peak_idx + 1:, cell] <= threshold)
        if not crossings.size:
            continue
        hi = peak_idx + 1 + int(crossings[0])
        lo = hi - 1
        v0, v1 = traces[lo, cell], traces[hi, cell]
        frac = 0.0 if v1 == v0 else (threshold-v0)/(v1-v0)
        crossing_time = times[lo] + frac*(times[hi]-times[lo])
        result[cell] = crossing_time - activation[cell]
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cpu_case", type=Path)
    parser.add_argument("cuda_case", type=Path)
    parser.add_argument("--write-interval", type=float, required=True)
    args = parser.parse_args()
    cpu_dirs, cuda_dirs = times_in(args.cpu_case), times_in(args.cuda_case)
    common_times = sorted(set(cpu_dirs) & set(cuda_dirs))
    if len(common_times) < 4:
        raise SystemExit("need at least four shared Vm output times")
    times = np.asarray(common_times, dtype=float)
    cpu_vm = np.stack([
        field(args.cpu_case, cpu_dirs[t], "Vm") for t in common_times
    ])
    cuda_vm = np.stack([
        field(args.cuda_case, cuda_dirs[t], "Vm") for t in common_times
    ])
    if cpu_vm.shape != cuda_vm.shape or not np.isfinite(cpu_vm).all() \
            or not np.isfinite(cuda_vm).all():
        raise SystemExit("Vm shape mismatch or nonfinite field values")
    cpu_activation = field(
        args.cpu_case, cpu_dirs[common_times[-1]], "activationTime", cpu_vm.shape[1]
    )
    cuda_activation = field(
        args.cuda_case, cuda_dirs[common_times[-1]], "activationTime", cuda_vm.shape[1]
    )
    cpu_active, cuda_active = cpu_activation > 0, cuda_activation > 0
    shared = cpu_active & cuda_active
    count_delta = abs(int(cpu_active.sum()) - int(cuda_active.sum()))
    mask_delta = int(np.count_nonzero(cpu_active ^ cuda_active))
    count_limit = max(
        ACTIVATED_COUNT_LIMIT_FLOOR,
        math.ceil(ACTIVATED_COUNT_LIMIT_FRACTION*cpu_vm.shape[1]),
    )
    act_diff = np.abs(cuda_activation[shared]-cpu_activation[shared])
    act_p95 = float(np.percentile(act_diff, 95))*1000 if act_diff.size else math.inf
    act_max = float(np.max(act_diff))*1000 if act_diff.size else math.inf

    shift = cuda_activation - cpu_activation
    cell_ids = np.arange(cpu_vm.shape[1])
    aligned_rows = []
    aligned_valid_rows = []
    for time_i, time in enumerate(times):
        query = time + shift
        valid = shared & (query >= times[0]) & (query <= times[-1])
        hi = np.clip(np.searchsorted(times, query, side="right"), 1, len(times)-1)
        lo = hi - 1
        span = times[hi] - times[lo]
        weight = np.divide(query-times[lo], span, out=np.zeros_like(query), where=span != 0)
        aligned_cuda = cuda_vm[lo, cell_ids]*(1-weight) + cuda_vm[hi, cell_ids]*weight
        aligned_rows.append(aligned_cuda - cpu_vm[time_i])
        aligned_valid_rows.append(valid)
    aligned_diff = np.asarray(aligned_rows)*1000.0
    aligned_valid = np.asarray(aligned_valid_rows)
    outside = np.abs(times[:, None]-cpu_activation[None, :]) > UPSTROKE_HALF_WIDTH_S
    mask = aligned_valid & outside
    aligned_rmse = float(np.sqrt(np.mean(aligned_diff[mask]**2))) if np.any(mask) else math.inf
    raw_diff_mv = (cuda_vm-cpu_vm)*1000
    raw_rmse = float(np.sqrt(np.mean(raw_diff_mv**2)))

    cpu_peak = np.max(cpu_vm[:, shared], axis=0) if np.any(shared) else np.array([])
    cuda_peak = np.max(cuda_vm[:, shared], axis=0) if np.any(shared) else np.array([])
    peak_p95 = (
        float(np.percentile(np.abs(cuda_peak-cpu_peak), 95))*1000
        if cpu_peak.size else math.inf
    )
    cpu_apd = apd90(times, cpu_vm, cpu_activation, shared)
    cuda_apd = apd90(times, cuda_vm, cuda_activation, shared)
    valid_apd = np.isfinite(cpu_apd) & np.isfinite(cuda_apd)
    if np.any(valid_apd):
        apd_diff_ms = np.abs(cuda_apd[valid_apd]-cpu_apd[valid_apd])*1000
        apd_p95 = float(np.percentile(apd_diff_ms, 95))
        ref_apd_ms = float(np.percentile(cpu_apd[valid_apd], 50))*1000
        apd_limit = max(APD90_ABS_LIMIT_MS, APD90_REL_LIMIT*ref_apd_ms)
    else:
        apd_p95, apd_limit = math.nan, math.inf

    print(f"shared times={len(times)} range={times[0]:g}..{times[-1]:g}s "
          f"sample_interval={args.write_interval:g}s cells={cpu_vm.shape[1]}")
    print(f"activated_cpu={int(cpu_active.sum())} activated_cuda={int(cuda_active.sum())} "
          f"count_delta={count_delta} mask_delta_cells={mask_delta} "
          f"mask_delta_limit={count_limit} "
          f"activation_p95_ms={act_p95:.6e} activation_max_ms={act_max:.6e}")
    print(f"Vm_raw_RMSE_mV={raw_rmse:.6e} "
          f"Vm_aligned_outside_upstroke_RMSE_mV={aligned_rmse:.6e} "
          f"peak_shift_p95_mV={peak_p95:.6e}")
    print(f"APD90_cells={int(valid_apd.sum())} APD90_shift_p95_ms={apd_p95:.6e} "
          f"APD90_limit_ms={apd_limit:.6e}")

    checks = {
        "activation p95": act_p95 <= ACTIVATION_P95_LIMIT_MS,
        "activated-cell mask": mask_delta <= count_limit,
        "aligned Vm RMSE outside upstroke":
            aligned_rmse <= ALIGNED_OUTSIDE_UPSTROKE_RMSE_LIMIT_MV,
        "peak voltage shift": peak_p95 <= PEAK_P95_LIMIT_MV,
    }
    if args.write_interval <= 0.001000001:
        checks["APD90 shift"] = (
            valid_apd.sum() > 0 and apd_p95 <= apd_limit
        )
    print("Acceptance: " + ", ".join(
        f"{name}={'PASS' if okay else 'FAIL'}" for name, okay in checks.items()
    ))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
