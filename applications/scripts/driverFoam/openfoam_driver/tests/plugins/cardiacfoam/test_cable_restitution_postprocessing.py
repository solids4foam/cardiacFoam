from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path

import pytest

from openfoam_driver.core.runtime.models import CaseConfig
from openfoam_driver.core.specs.spatial_pacing import generate_spatial_stimulus_lists
from openfoam_driver.plugins.cardiacfoam.tutorials.cable_1d_restitution import (
    _apply_case,
    _build_cases,
)
from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo


def _module():
    assert monorepo_root is not None
    path = (
        monorepo_root
        / "tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV"
        / "setup/postProcessing_cableRestitution.py"
    )
    spec = importlib.util.spec_from_file_location("cable_restitution_post", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Probe stations of the real 20 mm cable, in metres.
PROBE_X = (0.002, 0.005, 0.010, 0.015, 0.018)
RESTING_V = -0.084
PLATEAU_V = 0.030
APD_S = 0.300
CABLE_CV_M_PER_S = 3.36

# Measurements taken from the retained dt = 1e-6 s Stewart reference, used so
# the boundary-sweep regression tests below exercise the real timings.
REFERENCE_REPOLARIZATION90_S = 4.304086260869566
FIRST_SPONTANEOUS_S = 5.203520
S1_TIMES = [0.0, 1.0, 2.0, 3.0, 4.0]


def _write_trace(path: Path, beats, t_end: float = 6.0, dt: float = 5.0e-4) -> Path:
    """Square-wave action potentials on a uniform grid.

    ``beats`` is a list of ``(activation_time_at_first_probe, propagates)``.
    A non-propagating beat excites only the first probe, which is how local
    capture followed by conduction block appears in the probe traces.
    """
    header = "".join(f"# Probe {i} ({x} 5e-05 5e-05)\n" for i, x in enumerate(PROBE_X))
    lines = []
    for step in range(int(t_end / dt) + 1):
        time = step * dt
        values = []
        for index, x in enumerate(PROBE_X):
            value = RESTING_V
            for activation, propagates in beats:
                if not propagates and index > 0:
                    continue
                arrival = activation + (
                    (x - PROBE_X[0]) / CABLE_CV_M_PER_S if propagates else 0.0
                )
                if arrival <= time < arrival + APD_S:
                    value = PLATEAU_V
                    break
            values.append(value)
        lines.append(f"{time:.6f} " + " ".join(f"{v:g}" for v in values))
    path.write_text(header + "\n".join(lines) + "\n", encoding="ascii")
    return path


def _protocol(s2_times, **extra) -> dict:
    protocol = {
        "n_s1": len(S1_TIMES),
        "n_s2": len(s2_times),
        "s1_stimulus_times_s": list(S1_TIMES),
        "s2_stimulus_times_s": list(s2_times),
    }
    protocol.update(extra)
    return protocol


def _s1_beats():
    return [(time + 0.001, True) for time in S1_TIMES]


def _outcome(tmp_path: Path, beats, s2_times, name="Vm", **extra) -> dict:
    module = _module()
    trace = _write_trace(tmp_path / name, beats)
    return module.build_summary(trace, _protocol(s2_times, **extra))


# --------------------------------------------------------------------------
# Activation detection
# --------------------------------------------------------------------------


@skip_without_monorepo
def test_activation_crossing_is_interpolated() -> None:
    module = _module()
    crossings = module.detect_activation_crossings(
        [0.0, 0.1, 0.2, 0.3],
        [-80.0, 20.0, -80.0, -80.0],
    )
    assert crossings == pytest.approx([0.08])


@skip_without_monorepo
def test_plateau_notch_is_not_a_second_activation() -> None:
    """A dip to -1 and back is a notch, not a beat: the trace never recovered."""
    module = _module()
    crossings = module.detect_activation_crossings(
        [0.0, 0.1, 0.2, 0.21, 0.3],
        [-80.0, 20.0, -1.0, 1.0, -80.0],
    )
    assert crossings == pytest.approx([0.08])


@skip_without_monorepo
def test_second_activation_accepted_after_full_repolarization() -> None:
    """The same 130 ms separation counts once the trace returns to rest."""
    module = _module()
    crossings = module.detect_activation_crossings(
        [0.0, 0.1, 0.15, 0.2, 0.21, 0.3],
        [-80.0, 20.0, -80.0, -80.0, 1.0, -80.0],
    )
    assert len(crossings) == 2
    assert crossings[0] == pytest.approx(0.08)


# --------------------------------------------------------------------------
# Repolarization and APD
# --------------------------------------------------------------------------


@skip_without_monorepo
def test_apd_uses_beat_specific_mdp_and_interpolated_repolarization() -> None:
    module = _module()
    events = module.calculate_apd_events(
        [0.0, 0.08, 0.10, 0.20, 0.30, 0.40],
        [-80.0, -80.0, 20.0, 0.0, -60.0, -80.0],
        [0.08],
    )
    assert events[0]["mdp_V"] == -80.0
    assert events[0]["repolarization50_time_s"] == pytest.approx(0.25)
    assert events[0]["apd50_s"] == pytest.approx(0.17)
    assert events[0]["repolarization70_time_s"] == pytest.approx(0.2833333333333333)
    assert events[0]["repolarization90_level_V"] == -70.0
    assert events[0]["repolarization90_time_s"] == pytest.approx(0.35)
    assert events[0]["apd90_s"] == pytest.approx(0.27)
    assert events[0]["repolarization90_status"] == "measured"


@skip_without_monorepo
def test_repolarization_truncated_by_next_activation_is_labelled() -> None:
    """A premature next beat must not silently null the APD it cut short."""
    module = _module()
    events = module.calculate_apd_events(
        [0.0, 0.08, 0.10, 0.20, 0.28, 0.30],
        [-80.0, -80.0, 20.0, 0.0, 20.0, 20.0],
        [0.08, 0.27],
    )
    assert events[0]["repolarization90_time_s"] is None
    assert events[0]["apd90_s"] is None
    assert events[0]["repolarization90_status"] == "truncated_by_next_activation"


# --------------------------------------------------------------------------
# Stimulus association
# --------------------------------------------------------------------------


@skip_without_monorepo
def test_protocol_association_does_not_call_last_activation_s2() -> None:
    module = _module()
    protocol = {"s1_stimulus_times_s": [0.0, 1.0], "s2_stimulus_times_s": [1.5]}
    responses, spontaneous = module.associate_stimuli([0.01, 1.01, 1.20], protocol, 0.015)
    assert [response["captured"] for response in responses] == [True, True, False]
    assert spontaneous == [1.20]


@skip_without_monorepo
def test_association_window_widens_with_probe_distance() -> None:
    module = _module()
    positions = [(x, 0.0, 0.0) for x in PROBE_X]
    windows = [module.association_window_s(positions, i) for i in range(len(PROBE_X))]
    assert windows[0] == pytest.approx(module.STIMULUS_LATENCY_ALLOWANCE_S)
    assert windows == sorted(windows)
    # Bounded by the measured automaticity collision this protocol must survive.
    assert max(windows) < FIRST_SPONTANEOUS_S - (REFERENCE_REPOLARIZATION90_S + 0.850)


# --------------------------------------------------------------------------
# Protocol outcomes
# --------------------------------------------------------------------------


@skip_without_monorepo
def test_captured_and_propagated(tmp_path: Path) -> None:
    t_s2 = 4.500
    summary = _outcome(tmp_path, _s1_beats() + [(t_s2 + 0.0011, True)], [t_s2])
    assert summary["protocol_outcome"] == "captured_and_propagated"
    assert summary["central_cv"]["cv_m_per_s"] == pytest.approx(CABLE_CV_M_PER_S, rel=0.05)


@skip_without_monorepo
def test_stimulus_no_capture(tmp_path: Path) -> None:
    summary = _outcome(tmp_path, _s1_beats(), [4.500])
    assert summary["protocol_outcome"] == "stimulus_no_capture"
    assert "central_cv" not in summary


@skip_without_monorepo
def test_local_capture_propagation_block(tmp_path: Path) -> None:
    t_s2 = 4.500
    summary = _outcome(tmp_path, _s1_beats() + [(t_s2 + 0.0011, False)], [t_s2])
    assert summary["protocol_outcome"] == "local_capture_propagation_block"
    assert summary["restitution_metrics"][0]["s2_captured"] is True
    assert summary["restitution_metrics"][-1]["s2_captured"] is False
    assert "central_cv" not in summary


@skip_without_monorepo
def test_automaticity_before_s2(tmp_path: Path) -> None:
    t_s2 = 4.500
    summary = _outcome(tmp_path, _s1_beats() + [(4.400, True), (t_s2 + 0.0011, True)], [t_s2])
    assert summary["protocol_outcome"] == "automaticity_before_s2"
    assert summary["spontaneous_before_s2"] is True


@skip_without_monorepo
def test_competing_spontaneous_wave_after_valid_s2(tmp_path: Path) -> None:
    """A valid S2 followed by an unforced beat is not a clean capture."""
    t_s2 = 4.500
    summary = _outcome(
        tmp_path, _s1_beats() + [(t_s2 + 0.0011, True), (t_s2 + 0.400, True)], [t_s2]
    )
    assert summary["protocol_outcome"] == "competing_spontaneous_wave"
    assert summary["spontaneous_after_s2"] is True
    # The S2 wavefront was still measured cleanly, so CV survives the flag.
    assert summary["central_cv"]["cv_m_per_s"] == pytest.approx(CABLE_CV_M_PER_S, rel=0.05)


@skip_without_monorepo
def test_late_s2_does_not_absorb_the_measured_automaticity_beat(tmp_path: Path) -> None:
    """Regression for the real 850 ms boundary case.

    t(S2) = 5.154086 s and the measured first spontaneous activation is at
    5.203520 s -- 49.4 ms later. A flat 50 ms association window bound that
    unforced beat to S2 and reported a fabricated capture with a CV computed
    from automaticity.
    """
    t_s2 = REFERENCE_REPOLARIZATION90_S + 0.850
    assert FIRST_SPONTANEOUS_S - t_s2 == pytest.approx(0.049434, abs=1e-6)
    summary = _outcome(
        tmp_path,
        _s1_beats() + [(FIRST_SPONTANEOUS_S, True)],
        [t_s2],
        requested_di90_s=0.850,
        reference_repolarization90_s=REFERENCE_REPOLARIZATION90_S,
    )
    assert summary["protocol_outcome"] == "stimulus_no_capture"
    assert summary["spontaneous_after_s2"] is True
    assert "central_cv" not in summary


# --------------------------------------------------------------------------
# Restitution metrics
# --------------------------------------------------------------------------


@skip_without_monorepo
def test_s2_summary_keeps_coupling_activation_and_di90_distinct(tmp_path: Path) -> None:
    t_s2 = 4.500
    summary = _outcome(
        tmp_path,
        _s1_beats() + [(t_s2 + 0.0011, True)],
        [t_s2],
        requested_di90_s=0.1959,
        reference_repolarization90_s=REFERENCE_REPOLARIZATION90_S,
    )
    metric = summary["restitution_metrics"][0]
    assert summary["protocol_outcome"] == "captured_and_propagated"
    assert metric["stimulus_coupling_interval_s"] == pytest.approx(0.500)
    assert metric["activation_interval_s"] == pytest.approx(0.5001, abs=1e-3)
    # Coupling interval, activation interval and measured DI90 are three
    # different quantities and must not collapse onto one another.
    assert metric["measured_di90_s"] != pytest.approx(metric["stimulus_coupling_interval_s"])
    assert metric["measured_di90_s"] != pytest.approx(metric["activation_interval_s"])
    assert metric["measured_di90_status"] == "measured"


@skip_without_monorepo
def test_s2_apd_is_extracted_for_the_apd_restitution_axis(tmp_path: Path) -> None:
    t_s2 = 4.500
    summary = _outcome(tmp_path, _s1_beats() + [(t_s2 + 0.0011, True)], [t_s2])
    metric = summary["restitution_metrics"][0]
    assert metric["s2_apd90_s"] == pytest.approx(APD_S, abs=2e-3)
    assert metric["s2_apd70_s"] is not None
    assert metric["s2_apd50_s"] is not None
    assert metric["s1_apd90_s"] == pytest.approx(APD_S, abs=2e-3)


@skip_without_monorepo
def test_reference_repolarization_residual_is_recorded(tmp_path: Path) -> None:
    """The branch must be shown to reproduce the reference it was scheduled from.

    The synthetic S1 activates at 4.001 s with a 300 ms plateau, so the
    reference-probe repolarization90 lands at ~4.301 s.
    """
    t_s2 = 4.500
    beats = _s1_beats() + [(t_s2 + 0.0011, True)]

    matched = _outcome(
        tmp_path, beats, [t_s2], name="Vm_ok", reference_repolarization90_s=4.301
    )["reference_repolarization90_check"]
    assert matched["measured_s"] == pytest.approx(4.301, abs=1e-3)
    assert matched["consistent"] is True

    drifted = _outcome(
        tmp_path, beats, [t_s2], name="Vm_drift", reference_repolarization90_s=4.290
    )["reference_repolarization90_check"]
    assert drifted["expected_s"] == 4.290
    assert drifted["residual_s"] == pytest.approx(drifted["measured_s"] - 4.290)
    assert abs(drifted["residual_s"]) > drifted["tolerance_s"]
    assert drifted["consistent"] is False


@skip_without_monorepo
def test_reference_probe_and_central_segment_are_recorded(tmp_path: Path) -> None:
    t_s2 = 4.500
    summary = _outcome(tmp_path, _s1_beats() + [(t_s2 + 0.0011, True)], [t_s2])
    settings = summary["measurement_settings"]
    assert settings["reference_probe_index"] == 0
    assert settings["central_segment_probes"] == [1, 3]
    assert summary["central_cv"]["start_probe"] == 1
    assert summary["central_cv"]["end_probe"] == 3
    assert summary["restitution_metrics"][0]["is_reference_probe"] is True


# --------------------------------------------------------------------------
# Failure handling
# --------------------------------------------------------------------------


@skip_without_monorepo
def test_unreadable_probe_data_yields_analysis_failure_and_still_writes_rows(
    tmp_path: Path,
) -> None:
    module = _module()
    probes = tmp_path / "postProcessing" / "cableProbes" / "0"
    probes.mkdir(parents=True)
    (probes / "Vm").write_text("# Probe 0 (0.002 0 0)\n", encoding="ascii")
    output = tmp_path / "out"
    summary = module.write_summary_files(
        case_dir=tmp_path, output_dir=output, case_id="CASE_X"
    )
    assert summary["protocol_outcome"] == "analysis_failure"
    assert "No numeric probe data" in summary["analysis_failure_reason"]
    # A failed case must still leave a row behind, or it vanishes from the sweep.
    assert (output / "CASE_X_event_summary.json").is_file()
    assert (output / "CASE_X_events.csv").is_file()
    rows = (output / "CASE_X_restitution.csv").read_text().strip().splitlines()
    assert len(rows) == 2
    assert rows[1].startswith("CASE_X,analysis_failure")


@skip_without_monorepo
def test_censored_case_is_preserved_in_the_restitution_csv(tmp_path: Path) -> None:
    module = _module()
    case_dir = tmp_path / "case"
    probes = case_dir / "postProcessing" / "cableProbes" / "0"
    probes.mkdir(parents=True)
    t_s2 = 4.500
    _write_trace(probes / "Vm", _s1_beats() + [(4.400, True), (t_s2 + 0.0011, True)])
    (case_dir / ".cardiacfoam_protocol.json").write_text(
        json.dumps(_protocol([t_s2], requested_di90_s=0.5)), encoding="ascii"
    )
    output = tmp_path / "out"
    summary = module.write_summary_files(
        case_dir=case_dir, output_dir=output, case_id="CASE_C"
    )
    assert summary["protocol_outcome"] == "automaticity_before_s2"
    rows = (output / "CASE_C_restitution.csv").read_text().strip().splitlines()
    assert len(rows) == 1 + len(PROBE_X)
    assert all(row.startswith("CASE_C,automaticity_before_s2") for row in rows[1:])


# --------------------------------------------------------------------------
# No-S2 automaticity control
# --------------------------------------------------------------------------


@skip_without_monorepo
def test_no_s2_protocol_classifies_unassigned_activation_as_spontaneous(
    tmp_path: Path,
) -> None:
    module = _module()
    trace = _write_trace(tmp_path / "Vm", _s1_beats() + [(FIRST_SPONTANEOUS_S, True)])
    summary = module.build_summary(trace, _protocol([]))
    assert summary["protocol_outcome"] == "automaticity_observed"
    assert summary["spontaneous_activations_s"][0][0] == pytest.approx(
        FIRST_SPONTANEOUS_S, abs=1e-3
    )
    assert summary["automaticity_cycle_s"][0] == pytest.approx(
        FIRST_SPONTANEOUS_S - 4.001, abs=2e-3
    )
    assert "central_cv" not in summary


# --------------------------------------------------------------------------
# Driver-side scheduling
# --------------------------------------------------------------------------


def test_absolute_stimulus_schedule_preserves_fine_time() -> None:
    """%.6g would round this to 4.63409 s -- ~4 steps at deltaT = 1e-6 s."""
    schedule = generate_spatial_stimulus_lists(
        [0.0, 4.633857422617], "(0 0 0)", "(1 1 1)", "4e-3", "50000"
    )
    assert schedule["stimulusStartTimeList"] == "(0 4.63385742262)"


def _common_axes() -> dict:
    return {
        "ionic_models": ["Stewart"],
        "ionic_model_tissue_map": {"Stewart": ["myocyte"]},
        "dt_values": [0.001],
        "dx_values": [0.1],
        "solvers": ["implicit"],
        "conductivity_values": ["sigma"],
        "s2_intervals_ms": [500.0],
    }


def test_requested_di90_case_requires_and_records_repolarization_reference() -> None:
    with pytest.raises(ValueError, match="reference_repolarization90_s"):
        _build_cases(
            **_common_axes(),
            requested_di90_values_ms=[330.0],
            reference_repolarization90_s=None,
        )
    case = _build_cases(
        **_common_axes(),
        requested_di90_values_ms=[330.0],
        reference_repolarization90_s=REFERENCE_REPOLARIZATION90_S,
    )[0]
    assert case.case_id.endswith("_RDI90330")
    assert case.params["pacingMode"] == "requested_di90"
    assert case.params["requestedDI90_ms"] == 330.0
    # s2Interval means a coupling interval; a DI90 must not borrow that name.
    assert "s2Interval" not in case.params


def test_case_without_s2_is_not_named_after_a_pacing_value() -> None:
    case = _build_cases(**_common_axes(), n_s2=0)[0]
    assert case.case_id.endswith("_NOS2")


@skip_without_monorepo
def test_requested_di90_is_scheduled_as_an_absolute_stimulus_time(
    tmp_path: Path,
) -> None:
    """End-to-end: t(S2) = t(repolarization90_S1) + requestedDI90."""
    assert monorepo_root is not None
    case_root = tmp_path / "case"
    shutil.copytree(
        monorepo_root
        / "tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV",
        case_root,
    )
    case = _build_cases(
        **_common_axes(),
        requested_di90_values_ms=[330.0],
        reference_repolarization90_s=REFERENCE_REPOLARIZATION90_S,
    )[0]
    _apply_case(case_root, case, n_s1=5, n_s2=1)

    protocol = json.loads((case_root / ".cardiacfoam_protocol.json").read_text())
    expected_s2 = REFERENCE_REPOLARIZATION90_S + 0.330
    assert protocol["s2_stimulus_times_s"] == pytest.approx([expected_s2])
    assert protocol["s1_stimulus_times_s"] == pytest.approx(S1_TIMES)
    assert protocol["requested_di90_s"] == pytest.approx(0.330)
    assert protocol["reference_repolarization90_s"] == REFERENCE_REPOLARIZATION90_S
    # The coupling interval is a derived consequence here, not the input.
    assert protocol["s2_coupling_interval_s"] == pytest.approx(expected_s2 - 4.0)

    # Exactly five applied S1 stimuli and one S2 reach the case dictionary.
    written = (case_root / "constant" / "electroProperties").read_text()
    times = written.split("stimulusStartTimeList")[1].split("(")[1].split(")")[0].split()
    assert len(times) == 6
    assert float(times[-1]) == pytest.approx(expected_s2, abs=1e-9)
