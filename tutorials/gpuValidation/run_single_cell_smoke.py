#!/usr/bin/env python3
"""Run matched scalar and batched single-cell cases and compare their traces.

This is a diagnostic comparison. It measures voltage, APD90, and exported
state/current differences when those traces are available. It does not by
itself establish reference convergence, GPU correctness, or slab accuracy.
Each case uses a new output directory so existing results are preserved.
"""

import argparse
import csv
import shutil
import subprocess
from pathlib import Path

import numpy as np


MODELS = {
    "AlievPanfilov": ("myocyte", 0.4),
    "BuenoOrovio": ("epicardialCells", 2.0),
    "Courtemanche": ("myocyte", 60.0),
    "Fabbri": ("myocyte", 0.0),
    "Gaur": ("myocyte", 60.0),
    "Grandi": ("myocyte", 60.0),
    "PerisYague": ("myocyte", 60.0),
    "Stewart": ("myocyte", 60.0),
    "TNNP": ("epicardialCells", 60.0),
    "TWorld": ("epicardialCells", 60.0),
    "ToRORd_dynCl": ("epicardialCells", 60.0),
    "Trovato": ("myocyte", 60.0),
}


def control_dict(end_time, delta_t):
    return f"""FoamFile
{{
    version 2.0;
    format ascii;
    class dictionary;
    object controlDict;
}}
application cardiacFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {end_time};
deltaT {delta_t};
writeControl runTime;
writeInterval 100000;
writeFormat ascii;
writePrecision 9;
runTimeModifiable false;
"""


def electro_properties(model, tissue, amplitude, integrator, substeps, batched,
                       all_variables, write_frequency):
    model_name = f"{model}compactBatched" if batched else model
    batch_settings = (
        f"batchedIntegrator {integrator};\n    batchedSubsteps {substeps};"
        if batched else ""
    )
    export_spec = "()" if all_variables else "(Vm)"
    frequency_setting = (
        f"writeFrequency {write_frequency};" if write_frequency > 0 else ""
    )
    return f"""FoamFile
{{
    version 2.0;
    format ascii;
    class dictionary;
    object electroProperties;
}}
myocardiumSolver singleCellSolver;
singleCellSolverCoeffs
{{
    ionicModel {model_name};
    tissue {tissue};
    {batch_settings}
    solver RKF45;
    maxSteps 1000000000;
    singleCellStimulus
    {{
        stim_start 20;
        stim_duration 0.5;
        stim_amplitude {amplitude};
        stim_period_S1 1000;
        nstim1 3;
        stim_period_S2 0;
        nstim2 0;
    }}
    writeAfterTime 0;
    {frequency_setting}
    outputVariables
    {{
        ionic
        {{
            export {export_spec};
            debug ();
        }}
    }}
}}
"""


def run_case(case_dir, model, tissue, amplitude, args, batched, template):
    (case_dir / "constant").mkdir(parents=True)
    (case_dir / "system").mkdir()
    shutil.copy2(template / "constant/physicsProperties", case_dir / "constant")
    for name in ("blockMeshDict", "fvSchemes", "fvSolution"):
        shutil.copy2(template / "system" / name, case_dir / "system")
    (case_dir / "system/controlDict").write_text(
        control_dict(args.end_time, args.delta_t)
    )
    (case_dir / "constant/electroProperties").write_text(
        electro_properties(
            model, tissue, amplitude, args.integrator, args.substeps, batched,
            args.all_variables, args.write_frequency
        )
    )
    with (case_dir / "log.blockMesh").open("w") as log:
        mesh_result = subprocess.run(
            ["blockMesh", "-case", str(case_dir)], stdout=log, stderr=subprocess.STDOUT
        )
    if mesh_result.returncode:
        return {"status": f"blockMesh:{mesh_result.returncode}"}
    with (case_dir / "log.cardiacFoam").open("w") as log:
        result = subprocess.run(
            ["cardiacFoam", "-case", str(case_dir)],
            stdout=log, stderr=subprocess.STDOUT,
        )
    run_log = (case_dir / "log.cardiacFoam").read_text(errors="replace")
    backend = "scalar" if not batched else (
        "GPU" if "using CUDA device" in run_log else "CPU"
    )
    if result.returncode:
        return {"status": f"cardiacFoam:{result.returncode}", "backend": backend}
    if batched and args.require_gpu and backend != "GPU":
        return {"status": "GPU unavailable", "backend": backend}
    traces = list((case_dir / "postProcessing").glob("*.txt"))
    if len(traces) != 1:
        return {"status": f"expected one trace, found {len(traces)}", "backend": backend}
    trace = np.loadtxt(traces[0], skiprows=1)
    header = traces[0].read_text().splitlines()[0].split()
    if (trace.ndim != 2 or trace.shape[1] != len(header)
            or "Vm" not in header or not np.isfinite(trace).all()):
        return {"status": "invalid trace", "backend": backend}
    return {"status": "completed", "backend": backend,
            "trace": trace, "header": header}


def upward_crossing(trace, threshold):
    time, voltage = trace[:, 0], trace[:, 1]
    indices = np.flatnonzero((voltage[:-1] < threshold) & (voltage[1:] >= threshold))
    if not len(indices):
        return None
    i = indices[0]
    fraction = (threshold - voltage[i]) / (voltage[i + 1] - voltage[i])
    return float(time[i] + fraction * (time[i + 1] - time[i]))


def apd90_ms(trace):
    time, voltage = trace[:, 0], trace[:, 1]
    baseline = float(np.median(voltage[time < min(0.01, time[-1] / 5)]))
    peak_i = int(np.argmax(voltage))
    activation = upward_crossing(trace, -30.0)
    if activation is None:
        return None
    level = baseline + 0.1 * (voltage[peak_i] - baseline)
    indices = np.flatnonzero(
        (np.arange(len(voltage) - 1) >= peak_i)
        & (voltage[:-1] >= level)
        & (voltage[1:] < level)
    )
    if not len(indices):
        return None
    i = indices[0]
    fraction = (level - voltage[i]) / (voltage[i + 1] - voltage[i])
    repolarization = time[i] + fraction * (time[i + 1] - time[i])
    return float(1000 * (repolarization - activation))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--end-time", type=float, default=0.05)
    parser.add_argument("--delta-t", type=float, default=2e-6)
    parser.add_argument("--substeps", type=int, default=5)
    parser.add_argument("--integrator", choices=("euler", "rushLarsen"), default="rushLarsen")
    parser.add_argument("--all-variables", action="store_true")
    parser.add_argument("--write-frequency", type=float, default=0.0)
    parser.add_argument("--require-gpu", action="store_true")
    args = parser.parse_args()
    if (args.substeps < 1 or args.delta_t <= 0 or args.end_time <= 0
            or args.write_frequency < 0):
        parser.error("substeps, delta-t, and end-time must be positive")
    if args.all_variables and args.write_frequency == 0:
        args.write_frequency = 1e-4
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        parser.error(f"output directory already exists: {output_dir}")
    template = Path(__file__).resolve().parents[1] / "electrophysiologyProtocols/singleCell"
    output_dir.mkdir(parents=True)
    rows = []
    variable_rows = []
    for model in args.models:
        tissue, amplitude = MODELS[model]
        scalar = run_case(
            output_dir / "scalar" / model, model, tissue, amplitude, args, False, template
        )
        batched = run_case(
            output_dir / "batched" / model, model, tissue, amplitude, args, True, template
        )
        row = {
            "model": model,
            "tissue": tissue,
            "stim_amplitude": amplitude,
            "scalar_status": scalar["status"],
            "batched_status": batched["status"],
            "batched_backend": batched.get("backend", ""),
            "vm_rmse_mV": "",
            "vm_max_abs_mV": "",
            "peak_shift_ms": "",
            "activation_shift_ms": "",
            "vm_rmse_outside_upstroke_mV": "",
            "scalar_apd90_ms": "",
            "batched_apd90_ms": "",
            "apd90_difference_ms": "",
            "total_current": "",
            "integrated_current_error_over_abs_reference": "",
        }
        if scalar["status"] == batched["status"] == "completed":
            if scalar["header"] != batched["header"]:
                row["batched_status"] = "exported fields differ"
                rows.append(row)
                continue
            a_full, b_full = scalar["trace"], batched["trace"]
            vm_col = scalar["header"].index("Vm")
            a = a_full[:, [0, vm_col]]
            b = b_full[:, [0, vm_col]]
            if a.shape == b.shape and np.array_equal(a[:, 0], b[:, 0]):
                row["vm_rmse_mV"] = float(np.sqrt(np.mean((a[:, 1] - b[:, 1]) ** 2)))
                row["vm_max_abs_mV"] = float(np.max(np.abs(a[:, 1] - b[:, 1])))
                row["peak_shift_ms"] = float(
                    1000 * (b[np.argmax(b[:, 1]), 0] - a[np.argmax(a[:, 1]), 0])
                )
                activation_scalar = upward_crossing(a, -30.0)
                activation_batched = upward_crossing(b, -30.0)
                outside = np.ones(len(a), dtype=bool)
                if activation_scalar is not None:
                    outside &= np.abs(a[:, 0] - activation_scalar) > 0.002
                row["vm_rmse_outside_upstroke_mV"] = float(
                    np.sqrt(np.mean((a[outside, 1] - b[outside, 1]) ** 2))
                )
                if activation_scalar is not None and activation_batched is not None:
                    row["activation_shift_ms"] = float(
                        1000 * (activation_batched - activation_scalar)
                    )
                scalar_apd = apd90_ms(a)
                batched_apd = apd90_ms(b)
                if scalar_apd is not None:
                    row["scalar_apd90_ms"] = scalar_apd
                if batched_apd is not None:
                    row["batched_apd90_ms"] = batched_apd
                if scalar_apd is not None and batched_apd is not None:
                    row["apd90_difference_ms"] = batched_apd - scalar_apd
                if args.all_variables:
                    current_name = "Jion" if model == "BuenoOrovio" else "Iion_cm"
                    if current_name in scalar["header"]:
                        current_i = scalar["header"].index(current_name)
                        reference_current = a_full[:, current_i]
                        batched_current = b_full[:, current_i]
                        denominator = np.trapz(np.abs(reference_current), a_full[:, 0])
                        if denominator > 0:
                            row["total_current"] = current_name
                            row["integrated_current_error_over_abs_reference"] = float(
                                abs(np.trapz(batched_current - reference_current,
                                             a_full[:, 0])) / denominator
                            )
                    for i, field in enumerate(scalar["header"][1:], start=1):
                        error = a_full[:, i] - b_full[:, i]
                        reference = a_full[:, i]
                        scale = max(float(np.ptp(reference)),
                                    float(np.max(np.abs(reference))), 1e-12)
                        variable_rows.append({
                            "model": model,
                            "field": field,
                            "max_abs": float(np.max(np.abs(error))),
                            "rmse": float(np.sqrt(np.mean(error ** 2))),
                            "max_abs_over_scale": float(np.max(np.abs(error))) / scale,
                            "outside_upstroke_rmse": float(
                                np.sqrt(np.mean(error[outside] ** 2))
                            ),
                        })
            else:
                row["batched_status"] = "time grids differ"
        rows.append(row)
        print(f"{model}: scalar={row['scalar_status']} batched={row['batched_status']} "
              f"backend={row['batched_backend']}", flush=True)
    with (output_dir / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if args.all_variables and variable_rows:
        with (output_dir / "variable_summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(variable_rows[0]))
            writer.writeheader()
            writer.writerows(variable_rows)


if __name__ == "__main__":
    main()
