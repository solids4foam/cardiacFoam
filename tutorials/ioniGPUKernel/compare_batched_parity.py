#!/usr/bin/env python3
"""Compare CPU and CUDA batched slab fields, currents, and saved ODE states."""

import argparse
import re
import struct
from pathlib import Path

import numpy as np

from slab2D.setup.post_processing_slab_field_parity import field


STATE_SOURCE = {
    "TNNP": ("TNNP/TNNP_2004.H", "TNNP_STATES_NAMES"),
    "TWorld": ("TWorld/TWorld_2025.H", "TWorldSTATES_NAMES"),
}
MAGIC = 0x43524653
ABS_TOL = 1e-10


def state_names(model: str) -> list[str]:
    source, symbol = STATE_SOURCE[model]
    path = Path(__file__).resolve().parents[2] / "src/ionicModels" / source
    text = path.read_text()
    match = re.search(
        rf"static const char\*\s+{symbol}\s*\[[^]]+\]\s*=\s*\{{(.*?)\}};",
        text,
        re.DOTALL,
    )
    if match is None:
        raise ValueError(f"could not read state-name table {symbol} from {path}")
    return re.findall(r'"([^\"]+)"', match.group(1))


def read_restart_states(case: Path, time: str, model: str) -> np.ndarray:
    path = case / time / f"{model}compactBatchedState"
    data = path.read_bytes()
    if len(data) < 32:
        raise ValueError(f"short batched restart state file: {path}")
    magic, version, scalar_size, model_size, stride, rows = struct.unpack_from(
        "<IIIIQQ", data, 0
    )
    if magic != MAGIC or version != 1 or scalar_size != 8:
        raise ValueError(f"unsupported restart state header: {path}")
    offset = 32 + model_size
    file_model = data[32:offset].decode("ascii")
    if file_model != f"{model}compactBatched":
        raise ValueError(f"unexpected state model {file_model!r} in {path}")
    count = int(stride * rows)
    values = np.frombuffer(data, dtype="<f8", count=count, offset=offset)
    if values.size != count or not np.isfinite(values).all():
        raise ValueError(f"invalid or nonfinite restart state values: {path}")
    if offset + count * 8 != len(data):
        raise ValueError(f"unexpected trailing data in restart state file: {path}")
    return values.reshape((rows, stride))


def metrics(reference: np.ndarray, candidate: np.ndarray) -> tuple[float, float]:
    if reference.shape != candidate.shape:
        raise ValueError(f"field shapes differ: {reference.shape} != {candidate.shape}")
    difference = candidate - reference
    return float(np.sqrt(np.mean(difference**2))), float(np.max(np.abs(difference)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", choices=tuple(STATE_SOURCE))
    parser.add_argument("cpu_case", type=Path)
    parser.add_argument("gpu_case", type=Path)
    parser.add_argument("--times", nargs="+", default=["0.005", "0.01", "0.015"])
    args = parser.parse_args()
    names = state_names(args.model)
    passed = True

    for time in args.times:
        cpu_vm = field(args.cpu_case, time, "Vm")
        gpu_vm = field(args.gpu_case, time, "Vm")
        cpu_current = field(args.cpu_case, time, "ionicCurrent")
        gpu_current = field(args.gpu_case, time, "ionicCurrent")
        cpu_act = field(args.cpu_case, time, "activationTime", len(cpu_vm))
        gpu_act = field(args.gpu_case, time, "activationTime", len(gpu_vm))
        vm_rmse, vm_max = metrics(cpu_vm, gpu_vm)
        current_rmse, current_max = metrics(cpu_current, gpu_current)
        activation_diff = np.abs(cpu_act - gpu_act)
        cpu_active = cpu_act > 0
        gpu_active = gpu_act > 0
        shared = cpu_active & gpu_active
        p95 = (
            float(np.percentile(activation_diff[shared], 95))
            if np.any(shared)
            else float("nan")
        )
        cpu_states = read_restart_states(args.cpu_case, time, args.model)
        gpu_states = read_restart_states(args.gpu_case, time, args.model)
        if cpu_states.shape != gpu_states.shape:
            raise ValueError(f"restart state shapes differ at {time}s")
        if cpu_states.shape[1] != len(names):
            raise ValueError(
                f"state-name count {len(names)} does not match file stride "
                f"{cpu_states.shape[1]} for {args.model}"
            )

        print(
            f"{time}s cells={len(cpu_vm)} "
            f"Vm_RMSE={vm_rmse:.6e} Vm_max={vm_max:.6e} "
            f"IonicCurrent_RMSE={current_rmse:.6e} "
            f"IonicCurrent_max={current_max:.6e} "
            f"activation_cpu={int(cpu_active.sum())} "
            f"activation_gpu={int(gpu_active.sum())} "
            f"activation_p95_s={p95:.6e}"
        )

        if np.max(activation_diff) > ABS_TOL:
            passed = False
        if vm_max > ABS_TOL or current_max > ABS_TOL:
            passed = False
        for state_i, name in enumerate(names):
            state_rmse, state_max = metrics(
                cpu_states[:, state_i], gpu_states[:, state_i]
            )
            print(
                f"  {name}: RMSE={state_rmse:.6e} max={state_max:.6e}"
            )
            if state_max > ABS_TOL:
                passed = False

    print(
        f"CPU/GPU batched parity: {'PASS' if passed else 'FAIL'} "
        f"(absolute tolerance {ABS_TOL:g}; Vm/current/activation/state fields)"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
