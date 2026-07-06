from openfoam_driver.dashboard.plots import collect_plots


def test_collects_existing_pngs(tutorial_tree):
    case = tutorial_tree / "manufacturedSolutions" / "monodomainPseudoECG"
    (case / "postProcessing" / "ecg_plots.png").write_bytes(b"\x89PNG\r\n")
    found = collect_plots(case)
    assert any(p.name == "ecg_plots.png" for p in found)


def test_no_pngs_returns_empty(tutorial_tree):
    case = tutorial_tree / "PATHOS" / "RBBB"
    assert collect_plots(case) == []
