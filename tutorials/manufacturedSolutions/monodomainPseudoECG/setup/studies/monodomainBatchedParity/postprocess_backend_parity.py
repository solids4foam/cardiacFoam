#!/usr/bin/env python3
"""Compare two completed manufactured batched cases at matching output times."""

from __future__ import annotations

import argparse
import re
import struct
from pathlib import Path

import numpy as np


MAGIC = 0x43524653
MODEL_TYPE = "monodomainFDAManufacturedBatched"
STATE_NAMES = ("V", "u1", "u2", "u3")
ABS_TOL = 1.0e-10
FIELD_RE = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)", re.DOTALL
)


def field(case: Path, time_name: str, name: str, n_cells: int | None = None) -> np.ndarray:
    path = case / time_name / name
    data = path.read_text()
    match = FIELD_RE.search(data)
    if match is None:
        uniform = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", data)
        if uniform is None or n_cells is None:
            raise ValueError(f"expected scalar internal field: {path}")
        values = np.full(n_cells, float(uniform.group(1)))
    else:
        values = np.fromstring(match.group(2), sep=" ")
        if len(values) != int(match.group(1)):
            raise ValueError(f"invalid field length: {path}")
    if not np.isfinite(values).all():
        raise ValueError(f"nonfinite field: {path}")
    return values


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
    if values.size != count or offset + count * 8 != len(data):
        raise ValueError(f"incomplete state data: {path}")
    return values.reshape((rows, stride))


def compare_array(name: str, reference: np.ndarray, candidate: np.ndarray) -> bool:
    if reference.shape != candidate.shape:
        raise ValueError(f"{name} shapes differ: {reference.shape} != {candidate.shape}")
    difference = candidate - reference
    rmse = float(np.sqrt(np.mean(difference * difference)))
    maximum = float(np.max(np.abs(difference)))
    print(f"{name}: RMSE={rmse:.8e} max={maximum:.8e}")
    return maximum <= ABS_TOL


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_case", type=Path)
    parser.add_argument("candidate_case", type=Path)
    parser.add_argument("--times", nargs="+", default=["0.2"])
    args = parser.parse_args()
    passed = True

    for time_name in args.times:
        print(f"time={time_name}s")
        reference_vm = field(args.reference_case, time_name, "Vm")
        candidate_vm = field(args.candidate_case, time_name, "Vm")
        passed &= compare_array("Vm", reference_vm, candidate_vm)
        for name in ("ionicCurrent", "activationTime", "u1", "u2", "u3"):
            reference = field(args.reference_case, time_name, name, len(reference_vm))
            candidate = field(args.candidate_case, time_name, name, len(candidate_vm))
            passed &= compare_array(name, reference, candidate)

        reference_states = read_states(args.reference_case, time_name)
        candidate_states = read_states(args.candidate_case, time_name)
        for state_i, name in enumerate(STATE_NAMES):
            passed &= compare_array(
                f"state:{name}",
                reference_states[:, state_i],
                candidate_states[:, state_i],
            )

    print(
        f"Manufactured batched backend parity: {'PASS' if passed else 'FAIL'} "
        f"(max absolute tolerance {ABS_TOL:g})"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
