#!/usr/bin/env python3
"""Compare CPU and CUDA manufactured-batched final fields and ODE states."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np

from slab2D.setup.post_processing_slab_field_parity import field


MAGIC = 0x43524653
MODEL_TYPE = "monodomainFDAManufacturedBatched"
STATE_NAMES = ("V", "u1", "u2", "u3")
ABS_TOL = 1.0e-10


def read_states(case: Path, time_name: str) -> np.ndarray:
    path = case / time_name / f"{MODEL_TYPE}State"
    data = path.read_bytes()
    if len(data) < 32:
        raise ValueError(f"short restart state file: {path}")
    magic, version, scalar_size, model_size, stride, rows = struct.unpack_from(
        "<IIIIQQ", data, 0
    )
    if magic != MAGIC or version != 1 or scalar_size != 8:
        raise ValueError(f"unsupported restart state header: {path}")
    offset = 32 + model_size
    model_type = data[32:offset].decode("ascii")
    if model_type != MODEL_TYPE or stride != len(STATE_NAMES):
        raise ValueError(f"unexpected state header in {path}: {model_type}, {stride}")
    count = int(stride * rows)
    values = np.frombuffer(data, dtype="<f8", count=count, offset=offset)
    if values.size != count or not np.isfinite(values).all():
        raise ValueError(f"nonfinite or incomplete state data: {path}")
    if offset + count * 8 != len(data):
        raise ValueError(f"unexpected trailing bytes in {path}")
    return values.reshape((rows, stride))


def compare_array(name: str, cpu: np.ndarray, gpu: np.ndarray) -> bool:
    if cpu.shape != gpu.shape:
        raise ValueError(f"{name} shapes differ: {cpu.shape} != {gpu.shape}")
    difference = gpu - cpu
    rmse = float(np.sqrt(np.mean(difference * difference)))
    maximum = float(np.max(np.abs(difference)))
    print(f"{name}: RMSE={rmse:.8e} max={maximum:.8e}")
    return maximum <= ABS_TOL


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cpu_root", type=Path)
    parser.add_argument("gpu_root", type=Path)
    parser.add_argument("--scheme", choices=("godunov", "sbdf2"),
                        default="godunov")
    parser.add_argument("--times", nargs="+", default=["0.2"])
    args = parser.parse_args()
    passed = True

    for time_name in args.times:
        # All generated cases share the same mesh and final time. Select the
        # corresponding pair by their deltaT token in the directory name.
        cpu_cases = {p.name.split("_dt", 1)[1]: p
                     for p in args.cpu_root.glob(
                         f"cpu_{args.scheme}_*_dt*")}
        gpu_cases = {p.name.split("_dt", 1)[1]: p
                     for p in args.gpu_root.glob(
                         f"gpu_{args.scheme}_*_dt*")}
        if not cpu_cases or cpu_cases.keys() != gpu_cases.keys():
            raise ValueError("CPU/GPU manufactured case sets do not match")
        for dt_tag in sorted(cpu_cases, key=float, reverse=True):
            cpu_case = cpu_cases[dt_tag]
            gpu_case = gpu_cases[dt_tag]
            print(f"deltaT={dt_tag}, time={time_name}s")
            cpu_vm = field(cpu_case, time_name, "Vm")
            gpu_vm = field(gpu_case, time_name, "Vm")
            passed &= compare_array("Vm", cpu_vm, gpu_vm)
            for name in ("ionicCurrent", "activationTime", "u1", "u2", "u3"):
                cpu_values = field(cpu_case, time_name, name, len(cpu_vm))
                gpu_values = field(gpu_case, time_name, name, len(gpu_vm))
                passed &= compare_array(name, cpu_values, gpu_values)

            cpu_states = read_states(cpu_case, time_name)
            gpu_states = read_states(gpu_case, time_name)
            if cpu_states.shape != gpu_states.shape:
                raise ValueError("CPU/GPU restart state shapes differ")
            if cpu_states.shape[1] != len(STATE_NAMES):
                raise ValueError("manufactured state name count does not match data")
            for state_i, name in enumerate(STATE_NAMES):
                passed &= compare_array(
                    f"state:{name}", cpu_states[:, state_i],
                    gpu_states[:, state_i],
                )

    print(f"CPU/GPU manufactured-batched parity: {'PASS' if passed else 'FAIL'} "
          f"(max absolute tolerance {ABS_TOL:g})")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
