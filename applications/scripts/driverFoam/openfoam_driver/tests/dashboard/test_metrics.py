from pathlib import Path

from openfoam_driver.dashboard.metrics import load_dat, extract_metrics


def test_load_dat_parses_header_and_rows(tmp_path: Path):
    f = tmp_path / "pseudoECG.dat"
    f.write_text("# time  V1  V2\n0.0  1.0  2.0\n1.0  3.0  -4.0\n")
    cols, arr = load_dat(f)
    assert cols == ["time", "V1", "V2"]
    assert arr.shape == (2, 3)
    assert arr[1][2] == -4.0


def test_load_dat_ragged_rows_raise(tmp_path: Path):
    f = tmp_path / "bad.dat"
    f.write_text("# time V1\n0.0 1.0\n1.0\n")
    import pytest
    with pytest.raises(ValueError):
        load_dat(f)


def test_ecg_metrics_report_amplitude(tutorial_tree):
    case = tutorial_tree / "manufacturedSolutions" / "monodomainPseudoECG"
    m = extract_metrics("manufacturedSolutions/monodomainPseudoECG", case)
    # V1 spans 1..3 (pp=2), V2 spans -4..2 (pp=6) -> V2 is the largest lead
    assert m["n_samples"] == 2
    assert m["t_end"] == 1.0
    assert m["ecg_peak_to_peak"] == 6.0
    assert m["ecg_lead"] == "V2"


def test_unknown_case_returns_empty(tutorial_tree):
    case = tutorial_tree / "PATHOS" / "RBBB"
    assert extract_metrics("PATHOS/RBBB", case) == {}


def test_malformed_dat_yields_empty(tutorial_tree):
    case = tutorial_tree / "manufacturedSolutions" / "monodomainPseudoECG"
    (case / "postProcessing" / "pseudoECG.dat").write_text("# time V1\ngarbage\n")
    assert extract_metrics("manufacturedSolutions/monodomainPseudoECG", case) == {}
