from __future__ import annotations

def generate_spatial_s1_s2_stimulus_lists(
    s1_interval_ms: float, n_s1: int, s2_interval_ms: float, n_s2: int,
    bounds_min: str, bounds_max: str, duration_s: str, intensity: str
) -> dict[str, str]:
    times = []
    for i in range(n_s1):
        times.append(i * (s1_interval_ms / 1000.0))
    last_s1_time_s = (n_s1 - 1) * (s1_interval_ms / 1000.0)
    for i in range(n_s2):
        times.append(last_s1_time_s + (i + 1) * (s2_interval_ms / 1000.0))
    
    count = len(times)
    return {
        "stimulusStartTimeList": "(" + " ".join(f"{t:.6g}" for t in times) + ")",
        "stimulusLocationMinList": "(" + " ".join([bounds_min] * count) + ")",
        "stimulusLocationMaxList": "(" + " ".join([bounds_max] * count) + ")",
        "stimulusDurationList": "(" + " ".join([duration_s] * count) + ")",
        "stimulusIntensityList": "(" + " ".join([intensity] * count) + ")",
    }
