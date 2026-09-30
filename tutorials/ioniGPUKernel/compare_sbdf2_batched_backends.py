#!/usr/bin/env python3
"""Compare CPU-batched and CUDA-batched SBDF2 slab output and ODE states."""

import argparse
import struct
from pathlib import Path

import numpy as np

from slab2D.setup.post_processing_slab_field_parity import field


MAGIC = 0x43524653
VM_ABS_TOL = 1e-10
STATE_ABS_TOL = 1e-10
ACTIVATION_ABS_TOL_S = 1e-10
CURRENT_ABS_TOL = 1e-10
CURRENT_REL_TOL = 1e-12


def restart_states(case: Path, time: str, model: str) -> np.ndarray:
    path = case / time / f"{model}compactBatchedState"
    data = path.read_bytes()
    if len(data) < 32:
        raise ValueError(f"short compact-state file: {path}")
    magic, version, scalar_size, model_size, stride, rows = struct.unpack_from(
        "<IIIIQQ", data
    )
    offset = 32 + model_size
    state_model = data[32:offset].decode("ascii")
    if magic != MAGIC or version != 1 or scalar_size != 8:
        raise ValueError(f"unsupported compact-state header: {path}")
    if state_model != f"{model}compactBatched":
        raise ValueError(f"unexpected state model {state_model!r}: {path}")
    values = np.frombuffer(
        data, dtype="<f8", count=int(stride * rows), offset=offset
    )
    if values.size != stride * rows or not np.isfinite(values).all():
        raise ValueError(f"invalid or nonfinite compact states: {path}")
    if offset + values.size * 8 != len(data):
        raise ValueError(f"unexpected trailing bytes: {path}")
    return values.reshape((rows, stride))


def rmse_max(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    if a.shape != b.shape:
        raise ValueError(f"field dimensions differ: {a.shape} != {b.shape}")
    diff = np.abs(a - b)
    return float(np.sqrt(np.mean(diff * diff))), float(np.max(diff))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model")
    parser.add_argument("cpu_case", type=Path)
    parser.add_argument("gpu_case", type=Path)
    parser.add_argument("--times", nargs="+", default=["0.005", "0.01", "0.015"])
    parser.add_argument(
        "--scheme", choices=("godunov", "sbdf2"), default="sbdf2"
    )
    args = parser.parse_args()

    passed = True
    for time in args.times:
        cpu_vm = field(args.cpu_case, time, "Vm")
        gpu_vm = field(args.gpu_case, time, "Vm")
        cpu_current = field(args.cpu_case, time, "ionicCurrent")
        gpu_current = field(args.gpu_case, time, "ionicCurrent")
        cpu_activation = field(args.cpu_case, time, "activationTime", len(cpu_vm))
        gpu_activation = field(args.gpu_case, time, "activationTime", len(gpu_vm))
        vm_rmse, vm_max = rmse_max(cpu_vm, gpu_vm)
        current_rmse, current_max = rmse_max(cpu_current, gpu_current)
        activation_diff = np.abs(cpu_activation - gpu_activation)
        cpu_active = cpu_activation > 0
        gpu_active = gpu_activation > 0
        shared_active = cpu_active & gpu_active
        activation_p95 = (
            float(np.percentile(activation_diff[shared_active], 95))
            if np.any(shared_active)
            else float("nan")
        )
        cpu_states = restart_states(args.cpu_case, time, args.model)
        gpu_states = restart_states(args.gpu_case, time, args.model)
        state_rmse, state_max = rmse_max(cpu_states, gpu_states)
        current_scale = max(
            float(np.max(np.abs(cpu_current))),
            float(np.max(np.abs(gpu_current))),
        )
        current_tol = CURRENT_ABS_TOL + CURRENT_REL_TOL * current_scale

        print(
            f"{time}s cells={len(cpu_vm)} "
            f"Vm_RMSE={vm_rmse:.6e} Vm_max={vm_max:.6e} "
            f"Current_RMSE={current_rmse:.6e} Current_max={current_max:.6e} "
            f"Current_tol={current_tol:.6e} "
            f"activation_cpu={int(cpu_active.sum())} "
            f"activation_gpu={int(gpu_active.sum())} "
            f"activation_p95_s={activation_p95:.6e} "
            f"state_RMSE={state_rmse:.6e} state_max={state_max:.6e}"
        )
        passed &= vm_max <= VM_ABS_TOL
        passed &= current_max <= current_tol
        passed &= state_max <= STATE_ABS_TOL
        passed &= np.array_equal(cpu_active, gpu_active)
        passed &= float(np.max(activation_diff)) <= ACTIVATION_ABS_TOL_S

    print(
        f"CPU/GPU batched {args.scheme} parity: "
        f"{'PASS' if passed else 'FAIL'} "
        f"(Vm/state abs {VM_ABS_TOL:g}, current abs {CURRENT_ABS_TOL:g} "
        f"+ rel {CURRENT_REL_TOL:g}, activation abs {ACTIVATION_ABS_TOL_S:g} s)"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
