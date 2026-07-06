from openfoam_driver.dashboard.catalog import build_catalog


def test_catalog_merges_scan_status_metrics(tutorial_tree):
    cards = {c.case_id: c for c in build_catalog(tutorial_tree)}
    pe = cards["manufacturedSolutions/monodomainPseudoECG"]
    assert pe.status == "ran"
    assert pe.last_run is not None
    assert pe.metrics["ecg_lead"] == "V2"
    rbbb = cards["PATHOS/RBBB"]
    assert rbbb.status == "not_run"
    assert rbbb.is_3d is True


def test_catalog_marks_artifacts(tutorial_tree):
    case = tutorial_tree / "PATHOS" / "RBBB"
    (case / "heart.vtu").write_text("x")
    cards = {c.case_id: c for c in build_catalog(tutorial_tree)}
    assert "vtu" in cards["PATHOS/RBBB"].artifacts
