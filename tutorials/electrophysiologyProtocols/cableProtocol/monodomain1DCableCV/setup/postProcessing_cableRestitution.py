#!/usr/bin/env python3

from __future__ import annotations

import argparse
import bisect
import csv
import json
import re
import traceback
from pathlib import Path

# Vm is the solver's membrane potential in volts, not millivolts. The 0 V
# activation threshold and the relative repolarization levels are unit
# agnostic; every other absolute level below is derived from the trace's own
# amplitude so a model written in millivolts classifies identically.
VOLTAGE_UNIT = "V"
ACTIVATION_THRESHOLD_V = 0.0

# A threshold crossing only counts as a new activation once the trace has
# fallen back to this fraction of its own amplitude above the trace minimum.
# This replaces a former flat 50 ms debounce, which was wide enough to swallow
# a spontaneous beat arriving within 50 ms of a stimulated one. The remaining
# debounce exists only to reject numerical chatter on a single upstroke.
ACTIVATION_RESET_FRACTION = 0.25
ACTIVATION_DEBOUNCE_S = 0.002

# Stimulus-to-activation association window, evaluated per probe as
#     latency allowance + (probe distance from the first probe) / CV floor
# The CV floor sits ~2x below the slowest conduction velocity this protocol
# has measured (2.03 m/s), so genuinely decremental capture is still
# associated with its stimulus while an unforced beat arriving tens of
# milliseconds later is not. A flat 50 ms window previously bound the measured
# automaticity beat (t = 5.203520 s) to an S2 stimulus 49.4 ms earlier.
#
# The latency allowance covers stimulus-to-upstroke delay, which lengthens for
# a premature beat near refractoriness -- the regime the boundary sweep exists
# to probe. Both constants are bounded above by that same 49.4 ms collision:
# at these values the widest probe window is 31 ms, a 1.6x margin.
STIMULUS_LATENCY_ALLOWANCE_S = 0.015
ASSOCIATION_CV_FLOOR_M_PER_S = 1.0

PEAK_SEARCH_WINDOW_S = 0.10

# The probe whose measured repolarization90 defines the DI90 axis, and the
# probe pair whose separation defines the reported conduction velocity. They
# are deliberately different: repolarization90 varies by ~10 ms along this
# cable, so the DI90 axis is a proximal quantity while CV is a central-segment
# quantity, and a summary that does not say so is ambiguous.
REFERENCE_PROBE_INDEX = 0
CENTRAL_SEGMENT_PROBES = (1, 3)

# Tolerance on reproducing the conditioning reference used to schedule S2.
# Ten probe samples at the 1e-5 s probe write interval.
REFERENCE_REPOLARIZATION_TOLERANCE_S = 1.0e-4

REPOLARIZATION_PERCENTS = (50, 70, 90)

PROTOCOL_METADATA = ".cardiacfoam_protocol.json"


def parse_probe_file(path: Path) -> tuple[list[tuple[float, float, float]], list[float], list[list[float]]]:
    if not path.is_file():
        raise FileNotFoundError(f"Probe file not found: {path}")
    positions: list[tuple[float, float, float]] = []
    rows: list[list[float]] = []
    pattern = re.compile(r"#\s*Probe\s+\d+\s+\(([^)]+)\)")
    with path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            match = pattern.match(line)
            if match:
                coords = tuple(float(value) for value in match.group(1).split())
                if len(coords) != 3:
                    raise ValueError(f"Unexpected probe coordinates: {line}")
                positions.append(coords)
            elif line and not line.startswith("#"):
                rows.append([float(token) for token in line.split()])
    if not rows:
        raise ValueError(f"No numeric probe data found in {path}")
    times = [row[0] for row in rows]
    vm_data = [[row[column] for row in rows] for column in range(1, len(rows[0]))]
    if positions and len(positions) != len(vm_data):
        raise ValueError(f"Probe position/value mismatch: {len(positions)} vs {len(vm_data)}")
    return positions, times, vm_data


def interpolate_crossing(t0: float, y0: float, t1: float, y1: float, level: float) -> float:
    if y1 == y0:
        return t1
    return t0 + (level - y0) * (t1 - t0) / (y1 - y0)


def activation_reset_level(values: list[float]) -> float:
    """Absolute level the trace must fall below before re-arming activation.

    Derived from the trace's own amplitude so the rule is independent of
    whether Vm is expressed in volts or millivolts. Clamped to the activation
    threshold so a trace that never recovers that far still re-arms on a
    plain sub-threshold excursion rather than detecting only one beat.
    """
    low = min(values)
    high = max(values)
    return min(low + ACTIVATION_RESET_FRACTION * (high - low), ACTIVATION_THRESHOLD_V)


def detect_activation_crossings(times: list[float], values: list[float]) -> list[float]:
    reset_level = activation_reset_level(values)
    crossings: list[float] = []
    armed = True
    for index in range(1, len(times)):
        if values[index] < reset_level:
            armed = True
        if values[index - 1] < ACTIVATION_THRESHOLD_V <= values[index]:
            if not armed:
                continue
            crossing = interpolate_crossing(
                times[index - 1], values[index - 1], times[index], values[index], ACTIVATION_THRESHOLD_V
            )
            if crossings and crossing - crossings[-1] <= ACTIVATION_DEBOUNCE_S:
                continue
            crossings.append(crossing)
            armed = False
    return crossings


def calculate_apd_events(times: list[float], values: list[float], activations: list[float]) -> list[dict]:
    events: list[dict] = []
    previous_activation = times[0]
    for beat, activation in enumerate(activations):
        has_next_activation = beat + 1 < len(activations)
        next_activation = activations[beat + 1] if has_next_activation else times[-1]
        baseline_lo = bisect.bisect_left(times, previous_activation)
        baseline_hi = bisect.bisect_right(times, activation)
        peak_lo = bisect.bisect_left(times, activation)
        peak_hi = bisect.bisect_right(times, min(next_activation, activation + PEAK_SEARCH_WINDOW_S))
        if baseline_lo >= baseline_hi or peak_lo >= peak_hi:
            previous_activation = activation
            continue
        mdp = min(values[baseline_lo:baseline_hi])
        peak = max(values[peak_lo:peak_hi])
        peak_index = peak_lo + max(
            range(peak_hi - peak_lo), key=lambda offset: values[peak_lo + offset]
        )
        event = {
            "beat_index": beat,
            "activation_time_s": activation,
            "mdp_V": mdp,
            "peak_V": peak,
        }
        for percent in REPOLARIZATION_PERCENTS:
            level = mdp + (1.0 - percent / 100.0) * (peak - mdp)
            repolarization = None
            status = "not_reached"
            for index in range(peak_index + 1, len(times)):
                if has_next_activation and times[index] >= next_activation:
                    status = "truncated_by_next_activation"
                    break
                if values[index - 1] > level >= values[index]:
                    repolarization = interpolate_crossing(
                        times[index - 1], values[index - 1], times[index], values[index], level
                    )
                    status = "measured"
                    break
            event[f"repolarization{percent}_level_{VOLTAGE_UNIT}"] = level
            event[f"repolarization{percent}_time_s"] = repolarization
            event[f"repolarization{percent}_status"] = status
            event[f"apd{percent}_s"] = None if repolarization is None else repolarization - activation
        events.append(event)
        previous_activation = activation
    return events


def association_window_s(positions: list[tuple[float, float, float]], probe_index: int) -> float:
    """Widest plausible stimulus-to-activation delay at this probe."""
    if not positions:
        return STIMULUS_LATENCY_ALLOWANCE_S
    distance = abs(positions[probe_index][0] - positions[0][0])
    return STIMULUS_LATENCY_ALLOWANCE_S + distance / ASSOCIATION_CV_FLOOR_M_PER_S


def associate_stimuli(
    activations: list[float], protocol: dict, window_s: float
) -> tuple[list[dict], list[float]]:
    labelled = [
        (f"S1_{index + 1}", time)
        for index, time in enumerate(protocol.get("s1_stimulus_times_s", []))
    ] + [
        (f"S2_{index + 1}", time)
        for index, time in enumerate(protocol.get("s2_stimulus_times_s", []))
    ]
    unused = set(range(len(activations)))
    responses: list[dict] = []
    for label, stimulus_time in labelled:
        candidates = [
            index for index in unused
            if stimulus_time <= activations[index] <= stimulus_time + window_s
        ]
        selected = min(candidates, key=lambda index: activations[index]) if candidates else None
        if selected is not None:
            unused.remove(selected)
        responses.append({
            "stimulus_label": label,
            "stimulus_time_s": stimulus_time,
            "association_window_s": window_s,
            "captured": selected is not None,
            "activation_time_s": None if selected is None else activations[selected],
            # How many activations fell inside the window, and how much room
            # was left between the chosen one and the window's far edge. More
            # than one candidate means the window could not tell a stimulated
            # beat from something else arriving in the same interval.
            "candidate_count": len(candidates),
            "association_margin_s": (
                None if selected is None
                else stimulus_time + window_s - activations[selected]
            ),
        })
    return responses, [activations[index] for index in sorted(unused)]


def compute_segment_cv(positions: list[tuple[float, float, float]], activation_times: list[float]) -> list[dict]:
    segments: list[dict] = []
    for index in range(len(activation_times) - 1):
        dx = positions[index + 1][0] - positions[index][0]
        dt = activation_times[index + 1] - activation_times[index]
        if dt <= 0:
            raise ValueError(f"Non-positive activation delay between probes {index} and {index + 1}: {dt}")
        segments.append({
            "start_probe": index, "end_probe": index + 1,
            "dx_m": dx, "dt_s": dt, "cv_m_per_s": dx / dt,
        })
    return segments


def load_protocol_metadata(case_dir: Path) -> dict | None:
    path = case_dir / PROTOCOL_METADATA
    return json.loads(path.read_text(encoding="ascii")) if path.is_file() else None


def event_for_activation(events: list[dict], activation_time: float | None) -> dict | None:
    if activation_time is None or not events:
        return None
    return min(events, key=lambda event: abs(float(event["activation_time_s"]) - activation_time))


def _measurement_settings(positions: list[tuple[float, float, float]]) -> dict:
    return {
        "voltage_unit": VOLTAGE_UNIT,
        "activation_threshold": ACTIVATION_THRESHOLD_V,
        "activation_reset_fraction": ACTIVATION_RESET_FRACTION,
        "activation_debounce_s": ACTIVATION_DEBOUNCE_S,
        "stimulus_latency_allowance_s": STIMULUS_LATENCY_ALLOWANCE_S,
        "association_cv_floor_m_per_s": ASSOCIATION_CV_FLOOR_M_PER_S,
        "association_window_s_by_probe": [
            association_window_s(positions, index) for index in range(len(positions))
        ],
        "reference_probe_index": REFERENCE_PROBE_INDEX,
        "central_segment_probes": list(CENTRAL_SEGMENT_PROBES),
        "reference_repolarization90_tolerance_s": REFERENCE_REPOLARIZATION_TOLERANCE_S,
    }


def _reference_repolarization_check(protocol: dict, restitution_metrics: list[dict]) -> dict | None:
    expected = protocol.get("reference_repolarization90_s")
    if expected is None:
        return None
    measured = next(
        (
            metric.get("s1_repolarization90_time_s")
            for metric in restitution_metrics
            if metric["probe_index"] == REFERENCE_PROBE_INDEX
        ),
        None,
    )
    residual = None if measured is None else float(measured) - float(expected)
    return {
        "expected_s": float(expected),
        "measured_s": None if measured is None else float(measured),
        "residual_s": residual,
        "tolerance_s": REFERENCE_REPOLARIZATION_TOLERANCE_S,
        "consistent": (
            None if residual is None else abs(residual) <= REFERENCE_REPOLARIZATION_TOLERANCE_S
        ),
    }


def _summarize_no_s2(summary: dict, protocol: dict, apd_events: list[list[dict]], positions: list) -> None:
    spontaneous = summary["spontaneous_activations_s"]
    observed = any(spontaneous_probe for spontaneous_probe in spontaneous)
    summary["protocol_outcome"] = "automaticity_observed" if observed else "no_spontaneous_activation_observed"
    final_s1_label = f"S1_{int(protocol.get('n_s1', len(protocol.get('s1_stimulus_times_s', []))))}"
    final_s1_metrics = []
    automaticity_cycles = []
    for probe_index, (responses, probe_events) in enumerate(zip(summary["stimulus_responses"], spontaneous)):
        response = next((item for item in responses if item["stimulus_label"] == final_s1_label), None)
        activation = None if response is None else response["activation_time_s"]
        event = event_for_activation(apd_events[probe_index], activation)
        final_s1_metrics.append({
            "probe_index": probe_index,
            "activation_time_s": activation,
            "repolarization90_time_s": None if event is None else event["repolarization90_time_s"],
            "repolarization90_status": None if event is None else event["repolarization90_status"],
            "apd90_s": None if event is None else event["apd90_s"],
        })
        automaticity_cycles.append(
            None if activation is None or not probe_events else probe_events[0] - float(activation)
        )
    summary["final_s1_metrics"] = final_s1_metrics
    summary["automaticity_cycle_s"] = automaticity_cycles
    spontaneous_first = [events[0] for events in spontaneous if events]
    final_s1_activations = [
        metric["activation_time_s"] for metric in final_s1_metrics
        if metric["activation_time_s"] is not None
    ]
    if len(spontaneous_first) == len(positions) and len(final_s1_activations) == len(positions):
        probe_span = positions[-1][0] - positions[0][0]
        spontaneous_dt = spontaneous_first[-1] - spontaneous_first[0]
        paced_dt = float(final_s1_activations[-1]) - float(final_s1_activations[0])
        spontaneous_apparent_cv = None if spontaneous_dt <= 0 else probe_span / spontaneous_dt
        final_s1_cv = None if paced_dt <= 0 else probe_span / paced_dt
        summary["automaticity_spatial_metrics"] = {
            "probe_span_m": probe_span,
            "first_spontaneous_crossing_spread_s": max(spontaneous_first) - min(spontaneous_first),
            "first_to_last_probe_delay_s": spontaneous_dt,
            "apparent_first_to_last_probe_cv_m_per_s": spontaneous_apparent_cv,
            "final_s1_first_to_last_probe_cv_m_per_s": final_s1_cv,
            "apparent_to_final_s1_cv_ratio": (
                None if spontaneous_apparent_cv is None or final_s1_cv is None
                else spontaneous_apparent_cv / final_s1_cv
            ),
        }


def _summarize_s2(summary: dict, protocol: dict, apd_events: list[list[dict]]) -> list[float] | None:
    """Classify an S1-S2 branch and return the S2 activation times for CV.

    Outcome precedence, most specific first:
      analysis_failure           -- raised later, by the CV block
      automaticity_before_s2     -- an unforced beat pre-empted the S2 stimulus
      stimulus_no_capture        -- no probe responded to S2
      local_capture_propagation_block -- some but not all probes responded
      competing_spontaneous_wave -- S2 propagated, but an unforced beat followed
      captured_and_propagated    -- S2 propagated cleanly
    The two spontaneous_* booleans stay independent of the status so a
    partial-block case that also shows an unforced beat does not lose either
    fact to the single status field.
    """
    first_s2 = float(protocol["s2_stimulus_times_s"][0])
    spontaneous_before = any(
        any(time < first_s2 for time in probe) for probe in summary["spontaneous_activations_s"]
    )
    spontaneous_after = any(
        any(time >= first_s2 for time in probe) for probe in summary["spontaneous_activations_s"]
    )
    summary["spontaneous_before_s2"] = spontaneous_before
    summary["spontaneous_after_s2"] = spontaneous_after

    s2_responses = [
        next((item for item in probe_responses if item["stimulus_label"] == "S2_1"), None)
        for probe_responses in summary["stimulus_responses"]
    ]
    captured = [response is not None and response["captured"] for response in s2_responses]
    final_s1_label = f"S1_{int(protocol.get('n_s1', len(protocol.get('s1_stimulus_times_s', []))))}"

    restitution_metrics = []
    for probe_index, (probe_responses, s2_response) in enumerate(zip(summary["stimulus_responses"], s2_responses)):
        s1_response = next((item for item in probe_responses if item["stimulus_label"] == final_s1_label), None)
        s1_activation = None if s1_response is None else s1_response["activation_time_s"]
        s2_activation = None if s2_response is None else s2_response["activation_time_s"]
        s1_event = event_for_activation(apd_events[probe_index], s1_activation)
        s2_event = event_for_activation(apd_events[probe_index], s2_activation)
        repolarization90 = None if s1_event is None else s1_event["repolarization90_time_s"]
        if s2_activation is None:
            di90_status = "s2_not_captured"
        elif repolarization90 is not None:
            di90_status = "measured"
        elif s1_event is None:
            di90_status = "s1_event_missing"
        else:
            di90_status = s1_event["repolarization90_status"]
        metric = {
            "probe_index": probe_index,
            "is_reference_probe": probe_index == REFERENCE_PROBE_INDEX,
            "requested_di90_s": protocol.get("requested_di90_s"),
            "stimulus_coupling_interval_s": first_s2 - float(protocol["s1_stimulus_times_s"][-1]),
            "s2_captured": captured[probe_index],
            "activation_interval_s": (
                None if s1_activation is None or s2_activation is None
                else float(s2_activation) - float(s1_activation)
            ),
            "measured_di90_s": (
                None if repolarization90 is None or s2_activation is None
                else float(s2_activation) - float(repolarization90)
            ),
            "measured_di90_status": di90_status,
            "s1_activation_time_s": s1_activation,
            "s1_repolarization90_time_s": repolarization90,
            "s1_repolarization90_status": None if s1_event is None else s1_event["repolarization90_status"],
            "s1_apd90_s": None if s1_event is None else s1_event["apd90_s"],
            "s2_activation_time_s": s2_activation,
            "s2_stimulus_to_activation_latency_s": (
                None if s2_activation is None else float(s2_activation) - first_s2
            ),
        }
        for percent in REPOLARIZATION_PERCENTS:
            metric[f"s2_repolarization{percent}_time_s"] = (
                None if s2_event is None else s2_event[f"repolarization{percent}_time_s"]
            )
            metric[f"s2_repolarization{percent}_status"] = (
                None if s2_event is None else s2_event[f"repolarization{percent}_status"]
            )
            metric[f"s2_apd{percent}_s"] = None if s2_event is None else s2_event[f"apd{percent}_s"]
        restitution_metrics.append(metric)

    summary["restitution_metrics"] = restitution_metrics
    check = _reference_repolarization_check(protocol, restitution_metrics)
    if check is not None:
        summary["reference_repolarization90_check"] = check

    if spontaneous_before:
        summary["protocol_outcome"] = "automaticity_before_s2"
        return None
    if not any(captured):
        summary["protocol_outcome"] = "stimulus_no_capture"
        return None
    if not all(captured):
        summary["protocol_outcome"] = "local_capture_propagation_block"
        return None
    summary["protocol_outcome"] = (
        "competing_spontaneous_wave" if spontaneous_after else "captured_and_propagated"
    )
    return [float(response["activation_time_s"]) for response in s2_responses]


def build_summary(probe_path: Path, protocol: dict | None = None) -> dict:
    positions, times, vm_data = parse_probe_file(probe_path)
    activations = [detect_activation_crossings(times, values) for values in vm_data]
    apd_events = [calculate_apd_events(times, values, crossings) for values, crossings in zip(vm_data, activations)]
    summary: dict = {
        "probe_positions_m": [
            {"probe_index": i, "x": xyz[0], "y": xyz[1], "z": xyz[2]}
            for i, xyz in enumerate(positions)
        ],
        "measurement_settings": _measurement_settings(positions),
        "activation_events_s": activations,
        "apd_events": apd_events,
    }
    selected: list[float] | None = None
    if protocol is None:
        summary["selection_mode"] = "legacy_last_activation_no_protocol_metadata"
        if activations and all(crossings for crossings in activations):
            selected = [crossings[-1] for crossings in activations]
    else:
        summary["selection_mode"] = "explicit_protocol_stimulus_association"
        summary["protocol"] = protocol
        associations = [
            associate_stimuli(crossings, protocol, association_window_s(positions, probe_index))
            for probe_index, crossings in enumerate(activations)
        ]
        summary["stimulus_responses"] = [item[0] for item in associations]
        summary["spontaneous_activations_s"] = [item[1] for item in associations]
        summary["stimulus_association_ambiguous"] = any(
            response["candidate_count"] > 1
            for probe_responses in summary["stimulus_responses"]
            for response in probe_responses
        )
        if int(protocol.get("n_s2", 0)) == 0:
            _summarize_no_s2(summary, protocol, apd_events, positions)
        else:
            selected = _summarize_s2(summary, protocol, apd_events)

    if selected is not None and len(selected) == len(positions):
        required_probe = max(CENTRAL_SEGMENT_PROBES)
        try:
            if len(positions) <= required_probe:
                raise ValueError(
                    f"Central CV needs probes {CENTRAL_SEGMENT_PROBES}, "
                    f"but only {len(positions)} probes are present"
                )
            segments = compute_segment_cv(positions, selected)
            start, end = CENTRAL_SEGMENT_PROBES
            central_dx = positions[end][0] - positions[start][0]
            central_dt = selected[end] - selected[start]
            if central_dt <= 0:
                raise ValueError(f"Non-positive central activation delay: {central_dt}")
            summary["activation_times_s"] = selected
            summary["segments"] = segments
            summary["central_cv"] = {
                "start_probe": start, "end_probe": end, "dx_m": central_dx,
                "dt_s": central_dt, "cv_m_per_s": central_dx / central_dt,
            }
        except Exception as exc:  # noqa: BLE001 - recorded, never raised
            summary["protocol_outcome"] = "analysis_failure"
            summary["analysis_failure_reason"] = f"{type(exc).__name__}: {exc}"
            summary["activation_times_s"] = selected
    return summary


def _write_events_csv(path: Path, summary: dict) -> None:
    columns = ["probe_index", "beat_index", "activation_time_s", "mdp_V", "peak_V"]
    for percent in REPOLARIZATION_PERCENTS:
        columns += [
            f"repolarization{percent}_time_s",
            f"repolarization{percent}_status",
            f"apd{percent}_s",
        ]
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for probe_index, events in enumerate(summary.get("apd_events", [])):
            for event in events:
                row = [probe_index, event["beat_index"], event["activation_time_s"],
                       event["mdp_V"], event["peak_V"]]
                for percent in REPOLARIZATION_PERCENTS:
                    row += [
                        event[f"repolarization{percent}_time_s"],
                        event[f"repolarization{percent}_status"],
                        event[f"apd{percent}_s"],
                    ]
                writer.writerow(row)


def _write_restitution_csv(path: Path, summary: dict, case_id: str | None) -> None:
    """Long-form per-probe restitution row, including censored and failed cases."""
    metrics = summary.get("restitution_metrics") or []
    columns = [
        "case_id", "protocol_outcome", "probe_index", "is_reference_probe",
        "requested_di90_s", "stimulus_coupling_interval_s", "activation_interval_s",
        "measured_di90_s", "measured_di90_status", "s2_captured",
        "s1_activation_time_s", "s1_repolarization90_time_s", "s1_apd90_s",
        "s2_activation_time_s", "s2_stimulus_to_activation_latency_s",
        "s2_apd50_s", "s2_apd70_s", "s2_apd90_s", "s2_repolarization90_status",
        "central_cv_m_per_s",
    ]
    central_cv = (summary.get("central_cv") or {}).get("cv_m_per_s")
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        if not metrics:
            writer.writerow(
                [case_id, summary.get("protocol_outcome")] + [None] * (len(columns) - 2)
            )
            return
        for metric in metrics:
            writer.writerow(
                [case_id, summary.get("protocol_outcome")]
                + [metric.get(column) for column in columns[2:-1]]
                + [central_cv]
            )


def write_summary_files(*, case_dir: Path, output_dir: Path | None = None, case_id: str | None = None) -> dict:
    probe_files = list((case_dir / "postProcessing" / "cableProbes").glob("*/Vm"))
    if not probe_files:
        raise FileNotFoundError(f"No Vm probe files found in {case_dir}/postProcessing/cableProbes/")
    target_probe = max(probe_files, key=lambda path: float(path.parent.name))
    protocol = load_protocol_metadata(case_dir)
    try:
        summary = build_summary(target_probe, protocol)
    except Exception as exc:  # noqa: BLE001 - a failed case must still leave a row
        summary = {
            "protocol_outcome": "analysis_failure",
            "analysis_failure_reason": f"{type(exc).__name__}: {exc}",
            "analysis_failure_traceback": traceback.format_exc(),
            "probe_positions_m": [],
            "activation_events_s": [],
            "apd_events": [],
            "selection_mode": "analysis_failed_before_selection",
        }
        if protocol is not None:
            summary["protocol"] = protocol

    output_dir = output_dir or case_dir / "postProcessing"
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {"case_id": case_id, **summary}
    (output_dir / f"{case_id}_event_summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
    _write_events_csv(output_dir / f"{case_id}_events.csv", summary)
    _write_restitution_csv(output_dir / f"{case_id}_restitution.csv", summary, case_id)

    if "central_cv" in summary:
        cv_payload = {
            "case_id": case_id,
            "protocol_outcome": summary.get("protocol_outcome"),
            "selection_mode": summary.get("selection_mode"),
            "measurement_settings": summary.get("measurement_settings"),
            "probe_positions_m": summary.get("probe_positions_m"),
            "activation_times_s": summary.get("activation_times_s"),
            "segments": summary.get("segments"),
            "central_cv": summary["central_cv"],
        }
        (output_dir / f"{case_id}_cv_summary.json").write_text(json.dumps(cv_payload, indent=2) + "\n", encoding="ascii")
        with (output_dir / f"{case_id}_activation_times.csv").open("w", newline="", encoding="ascii") as handle:
            writer = csv.writer(handle)
            writer.writerow(["probe_index", "x_mm", "selected_activation_time_ms"])
            for index, activation in enumerate(summary["activation_times_s"]):
                writer.writerow([index, 1e3 * positions_x(summary, index), 1e3 * activation])
        with (output_dir / f"{case_id}_segments.csv").open("w", newline="", encoding="ascii") as handle:
            writer = csv.DictWriter(handle, fieldnames=["start_probe", "end_probe", "dx_m", "dt_s", "cv_m_per_s"])
            writer.writeheader()
            writer.writerows(summary["segments"])
    return summary


def positions_x(summary: dict, index: int) -> float:
    return float(summary["probe_positions_m"][index]["x"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", default=".")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--case-id", default=None)
    args = parser.parse_args()
    case_dir = Path(args.case_dir).resolve()
    case_id = args.case_id
    if case_id is None:
        sentinel = case_dir / ".driverfoam_case_id"
        if not sentinel.exists():
            raise SystemExit(f"--case-id not supplied and sentinel {sentinel} not found")
        case_id = sentinel.read_text(encoding="ascii").strip()
    write_summary_files(case_dir=case_dir, output_dir=Path(args.output_dir).resolve(), case_id=case_id)
